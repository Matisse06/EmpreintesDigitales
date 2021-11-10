package cs107;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Provides tools to compare fingerprint.
 */
public class Fingerprint {


	/**
	 * The number of pixels to consider in each direction when doing the linear
	 * regression to compute the orientation.
	 */
	public static final int ORIENTATION_DISTANCE = 16;

	/**
	 * The maximum distance between two minutiae to be considered matching.
	 */
	public static final int DISTANCE_THRESHOLD = 5;

	/**
	 * The number of matching minutiae needed for two fingerprints to be considered
	 * identical.
	 */
	public static final int FOUND_THRESHOLD = 20;

	/**
	 * The distance between two angle to be considered identical.
	 */
	public static final int ORIENTATION_THRESHOLD = 20;

	/**
	 * The offset in each direction for the rotation to test when doing the
	 * matching.
	 */
	public static final int MATCH_ANGLE_OFFSET = 2;

	/**
	 * Returns an array containing the value of the 8 neighbours of the pixel at
	 * coordinates <code>(row, col)</code>.
	 * <p>
	 * The pixels are returned such that their indices corresponds to the following
	 * diagram:<br>
	 * ------------- <br>
	 * | 7 | 0 | 1 | <br>
	 * ------------- <br>
	 * | 6 | _ | 2 | <br>
	 * ------------- <br>
	 * | 5 | 4 | 3 | <br>
	 * ------------- <br>
	 * <p>
	 * If a neighbours is out of bounds of the image, it is considered white.
	 * <p>
	 * If the <code>row</code> or the <code>col</code> is out of bounds of the
	 * image, the returned value should be <code>null</code>.
	 *
	 * @param image array containing each pixel's boolean value.
	 * @param row   the row of the pixel of interest, must be between
	 *              <code>0</code>(included) and
	 *              <code>image.length</code>(excluded).
	 * @param col   the column of the pixel of interest, must be between
	 *              <code>0</code>(included) and
	 *              <code>image[row].length</code>(excluded).
	 * @return An array containing each neighbours' value.
	 */

//determine la couleur du voisin (methode utilise par la suite dans getNeighbours
	public static boolean color(boolean[][] image, int row, int col) {
		boolean couleur = true;
		if( (row >= image.length || row < 0) || (col >= image[0].length || col < 0 ) ) {
			couleur = false;
		} else {
			couleur = image[row][col];
		}
		return couleur;
	}

	public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
		assert (image != null);
		if ((row < 0 || row >= image.length) || (col < 0 || col >= image[0].length)) {
			System.out.println("Pixel donne n'est pas dans l'image");
			return null;
		}

		boolean[] neighbours = new boolean[8];
		neighbours[0] = color(image, row-1, col);
		for(int i = -1; i < 2; ++i) {
			neighbours[i+2] = color(image, row+i, col+1);
		}
		neighbours[4] = color(image, row+1, col);
		neighbours[5] = color(image, row+1, col-1);
		neighbours[6] = color(image, row, col-1);
		neighbours[7] = color(image, row-1, col-1);

		return neighbours;



		// special case that is not expected (the image is supposed to have been checked
		// earlier)
		//TODO implement
	}

	/**
	 * Computes the number of black (<code>true</code>) pixels among the neighbours
	 * of a pixel.
	 *
	 * @param neighbours array containing each pixel value. The array must respect
	 *                   the convention described in
	 *                   {@link #getNeighbours(boolean[][], int, int)}.
	 * @return the number of black neighbours.
	 */

	public static int blackNeighbours(boolean[] neighbours) {

		assert (neighbours != null);

		int black = 0;

		for (int i = 0; i < neighbours.length; ++i) {
			if (neighbours[i]) {
				++black;
			}
		}

		return black;

	}



	/**

	 * Computes the number of white to black transitions among the neighbours of

	 * pixel.

	 *

	 * @param neighbours array containing each pixel value. The array must respect

	 *                   the convention described in

	 *                   {@link #getNeighbours(boolean[][], int, int)}.

	 * @return the number of white to black transitions.

	 */

	public static int transitions(boolean[] neighbours) {

		assert (neighbours != null);

		int numberOfTransitions = 0;

		for (int i = 0; i < neighbours.length; ++i) {
			int j = i + 1;

			if (i == 7) {
				j = 0;
			}
			if (!neighbours[i] && neighbours[j]) {
				++ numberOfTransitions;
			}
		}


		return numberOfTransitions;

	}


	/**

	 * Returns <code>true</code> if the images are identical and false otherwise.

	 *

	 * @param image1 array containing each pixel's boolean value.
	 * @param image2 array containing each pixel's boolean value.

	 * @return <code>True</code> if they are identical, <code>false</code>

	 *         otherwise.

	 */

	public static boolean identical(boolean[][] image1, boolean[][] image2) {

		assert ((image1 != null) && (image2 != null));

		//checks if two tableaux are identical
		if ((image1.length != image2.length) || (image1[0].length != image2[0].length) )
			return false;

		for (int i = 0; i < image1.length; ++i) {
			for (int j = 0; j < image1[i].length; ++j) {

				if ((image1[i][j]) != (image2[i][j])) {
					return false;
				}
			}
		}
		return true;
	}

	/**

	 * Internal method used by {@link #thin(boolean[][])}.
	 *
	 * @param //image array containing each pixel's boolean value.
	 * @param step  the step to apply, Step 0 or Step 1.

	 * @return A new array containing each pixel's value after the step.

	 */

	public static boolean[][] thinningStep(boolean[][] imageTest, int step) {

		assert ((step == 1 || step == 2) && (imageTest != null));

		boolean[][] imageTest4 = new boolean[imageTest.length][imageTest[0].length];

		for (int i = 0; i < imageTest.length; ++i) {
			for (int j = 0; j < imageTest[i].length; ++j) {
				imageTest4[i][j] = imageTest[i][j];
			}
		}

// boucle for pour déterminer si chaque pixel remplit les conditions

		for (int i = 0; i < imageTest.length; ++i) {
			for (int j = 0; j < imageTest[i].length; ++j) {

				// tests communs aux step 1 et 2

				boolean pixelNoir = imageTest[i][j];
				boolean[] neighbours = getNeighbours(imageTest, i, j);

				if ((pixelNoir) && (neighbours != null) &&
						(blackNeighbours(neighbours)) >= 2 && (blackNeighbours(neighbours) <= 6) &&
						(transitions(neighbours) == 1)) {

					// test step 1

					if (step == 0) {

						if (((getNeighbours(imageTest, i, j)[0] == false) ||(getNeighbours(imageTest, i, j)[2] == false) || (getNeighbours(imageTest, i, j)[4] == false)) &&
								((getNeighbours(imageTest, i, j)[2] == false) ||(getNeighbours(imageTest, i, j)[4] == false) || (getNeighbours(imageTest, i, j)[6] == false))) {
							imageTest4[i][j] = false;
						}

						// test step 2

					} else if (step == 1) {

						if (((getNeighbours(imageTest, i, j)[0] == false) ||(getNeighbours(imageTest, i, j)[2] == false) || (getNeighbours(imageTest, i, j)[6] == false)) &&
								((getNeighbours(imageTest, i, j)[0] == false) ||(getNeighbours(imageTest, i, j)[4] == false) || (getNeighbours(imageTest, i, j)[6] == false))) {
							imageTest4[i][j] = false;
						}
					}
				}
			}
		}
		return imageTest4;
	}

	/**

	 * Compute the skeleton of a boolean image.

	 * @param image array containing each pixel's boolean value.
	 * @return array containing the boolean value of each pixel of the image after
	 *         applying the thinning algorithm.

	 */

	public static boolean[][] thin(boolean[][] image) {

		assert (image != null);

		// création tableau test 1 pour ne pas modifier le tableau initial
		// création tableau test 2 qui va être modifié (step 1)
		// création tableau test 3 qui va être modifié (step 2)

		boolean[][] imageTest1 = new boolean[image.length][image[0].length];

		for (int i = 0; i < image.length; ++i) {
			for (int j = 0; j < image[i].length; ++j) {
				imageTest1[i][j] = image[i][j];
			}
		}

		boolean change;

		do {
			boolean[][] imageTest2 = thinningStep(imageTest1, 0);
			boolean[][] imageTest3 = thinningStep(imageTest2, 1);

			if (!identical(imageTest1, imageTest3)) {
				for (int i = 0; i < imageTest1.length; ++i) {
					for (int j = 0; j < imageTest1[i].length; ++j) {
						imageTest1[i][j] = imageTest3[i][j];
					}
				}
			}
			//

			change = !identical(imageTest2, imageTest3);


		} while (change);

		return imageTest1;

	}


	/**
	 * Compute the skeleton of a boolean image.
	 *
	 * @param //image array containing each pixel's boolean value.
	 * @return array containing the boolean value of each pixel of the image after
	 *         applying the thinning algorithm.
	 */

  /**
   * Computes all pixels that are connected to the pixel at coordinate
   * <code>(row, col)</code> and within the given distance of the pixel.
   *
   * @param //image    array containing each pixel's boolean value.
   * @param row      the first coordinate of the pixel of interest.
   * @param col      the second coordinate of the pixel of interest.
   * @param //distance the maximum distance at which a pixel is considered.
   * @return An array where <code>true</code> means that the pixel is within
   *         <code>distance</code> and connected to the pixel at
   *         <code>(row, col)</code>.
   */

	public static int rowNeighbour(int i, int row, int col) { // determine la ligne du pixel voisin de getNeighbours
		if(i==0 || i==1 || i==7) {
			return row-1;
		}
		if(i==2 || i==6) {
			return row;
		}
		if(i==3 || i==4 || i==5) {
			return row+1;
		}
		return -1;
	}

	public static int colNeighbour(int i, int row, int col) { // determine la colonne du pixel voisin de getNeighbours
		if(i==0 || i==4) {
			return col;
		}
		if(i==1 || i==2 || i==3) {
			return col+1;
		}
		if(i==5 || i==6 || i==7)  {
			return col-1;
		}
		return -1;
	}

	/**
	 *
	 * @param a
	 * @param b
	 * @return
	 */


	/**
	 *
	 * @param image
	 * @param row
	 * @param col
	 * @param distance
	 * @return
	 */
	public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {
		boolean[][] connected = new boolean[image.length][image[0].length];
		connected[row][col] = true;

		// idee 1 :
		boolean finit = false;
		while (finit == false) {
			finit = true;
			for (int i = 0; i < image.length; i++) {
				for (int j = 0; j < image[0].length; ++j) {
					if (image[i][j] == false) {
						//pixel blanc, go pixel suivant
					} else {
						if ((i < row - distance || i > row + distance) || (j < col - distance || j > col + distance)) {
							// trop loin de la minutiae, go suivant
						} else {
							if (connected[i][j] == true) {
								// si pixel noir deja dans tableau connected, go pixel suivant
							} else {
								boolean[] neighbours = getNeighbours(image, i, j);
								for (int g = 0; g < neighbours.length; ++g) {
									if (neighbours[g] == true) {
										int rowPixNoir = rowNeighbour(g, i, j);                // determiner position pixel noir voisin
										int colPixNoir = colNeighbour(g, i, j);
										if (connected[rowPixNoir][colPixNoir] == true) {        // si le pixel voisin est lie a minutiae good!
											connected[i][j] = true;
											finit = false;
											// on a trouve un pixel lie, recommencer processus
										}
									}
								}
							}
						}
					}

				}
			}

		}
		return connected;
	}

		/**
		 * Computes the slope of a minutia using linear regression.
		 *
		 * @param connectedPixels the result of
		 *                        {@link #connectedPixels(boolean[][], int, int, int)}.
		 * @param //row             the row of the minutia.
		 * @param //col             the col of the minutia.
		 * @return the slope.
		 */
		public static double computeSlope ( boolean[][] connectedPixels, int rowM, int colM){

			assert (connectedPixels != null);

			// initialisation des variables

			double slope = 0.;

			double somme_x2 = 0;
			double somme_y2 = 0;
			double somme_xy = 0;


			for (int row = 0; row < connectedPixels.length; ++row) {
				for (int col = 0; col < connectedPixels[row].length; ++col) {

					// x et y exprimés en fonction de l'origine du repère, centré sur la minutie
					int x = col - colM;
					int y = rowM - row;

					if (connectedPixels[row][col]) {
						somme_x2 += x * x;
						somme_y2 += y * y;
						somme_xy += x * y;
					}
				}
			}

			if (somme_x2 >= somme_y2) {
				slope = somme_xy / somme_x2;
			} else {
				slope = somme_y2 / somme_xy;
			}
		
			if (somme_x2 == 0) {
				slope = Double.POSITIVE_INFINITY;
				return slope;
			} 
			
			return slope;
		}


		/**
		 * Computes the orientation of a minutia in radians.
		 *
		 * @param connectedPixels the result of
		 *                        {@link #connectedPixels(boolean[][], int, int, int)}.
		 * @param //row             the row of the minutia.
		 * @param //col             the col of the minutia.
		 * @param slope           the slope as returned by
		 *                        {@link #computeSlope(boolean[][], int, int)}.
		 * @return the orientation of the minutia in radians.
		 */
		
		
		public static double computeAngle ( boolean[][] connectedPixels, int rowM, int colM, double slope){
			

			// calcul de l'angle (rad) formé par le vecteur direction
			double angle = Math.atan(slope);

			// variables qui comptent le nombre de pixels au dessus et en dessous
			int pixelUp = 0;
			int pixelDown = 0;

			// boucle qui vérifie que les conditions soient respectées
			
			for (int row = 0; row < connectedPixels.length; ++row) {
				for (int col = 0; col < connectedPixels[row].length; ++col) {

					// x et y exprimés en fonction de l'origine du repère, centré sur la minutie
					int x = col - colM;
					int y = rowM - row;

					boolean pixelNoir = connectedPixels[row][col];

					if (pixelNoir && y >= ((-1 / slope) * x)) {
						pixelUp += 1;
					} else if (pixelNoir) {
						pixelDown += 1;
					}
				}
			}
		
		
			// conditions qui déterminent l'ajout potentiel de PI à l'angle calculé au début
			if ((slope == Double.POSITIVE_INFINITY && (pixelUp > pixelDown))) {
				angle = Math.PI / 2;

			} else if ((slope == Double.POSITIVE_INFINITY && (pixelUp < pixelDown))) {
				angle = -Math.PI / 2;

			} else if ((angle > 0 && (pixelDown > pixelUp)) || (angle < 0 && (pixelDown <= pixelUp))) {
				angle += Math.PI;
			}
			
			

			// retourne l'angle en radian

			
			return angle;
		}

		/**
		 * Computes the orientation of the minutia that the coordinate <code>(row,
		 * col)</code>.
		 *
		 * @param image    array containing each pixel's boolean value.
		 * @param row      the first coordinate of the pixel of interest.
		 * @param col      the second coordinate of the pixel of interest.
		 * @param distance the distance to be considered in each direction to compute
		 *                 the orientation.
		 * @return The orientation in degrees.
		 */
		public static int computeOrientation ( boolean[][] image, int row, int col, int distance){

				// angle en radian
				boolean[][] neighbours = connectedPixels(image, row, col, ORIENTATION_DISTANCE);
				double slope = computeSlope(neighbours, row, col);
				double angleR = computeAngle(neighbours, row, col, slope);

				// conversion de l'angle en degré
				angleR = Math.toDegrees(angleR);

				// angle en degré (rounded)
				int angleD = (int) Math.round(angleR);

				if (angleD < 0) {
					angleD += 360;
				}
				
				return angleD;
			}


		/**
		 * Extracts the minutiae from a thinned image..
		 *
		 * @param image array containing each pixel's boolean value.
		 * @return The list of all minutiae. A minutia is represented by an array where
		 *         the first element is the row, the second is column, and the third is
		 *         the angle in degrees.
		 * @see #thin(boolean[][])
		 *
		 */
		public static List<int[]> extract(boolean[][] image){
			ArrayList<int[]> listeMinutie = new ArrayList<int[]>();   //creation of list with minutiaes.

			for(int i = 1 ; i < image.length - 1 ; ++i){
				for(int j = 1; j < image[0].length - 1 ; ++j) {     	// scan of image, starting at pixel 1 and ending at col-1
					if(image[i][j]) {								// if black pixel (because minutiaes are black)
						int transi = transitions(getNeighbours(image, i, j));    //calculating transitions
						if(transi == 1 || transi == 3) {			// if terminaison or bifurcation type of minutiae detected
							int[] minutie = new int[3]; // tab with all informations for 1 minutiae
							minutie[0] = i;
							minutie[1] = j;
							minutie[2] = computeOrientation(image, i, j, ORIENTATION_DISTANCE);
							listeMinutie.add(minutie);
						}
					}
				}
			}
			return listeMinutie;	// returns list of minuatiaes( their location + orientation)
		}

		/**
		 * Applies the specified rotation to the minutia.
		 *
		 * @param minutia   the original minutia.
		 * @param centerRow the row of the center of rotation.
		 * @param centerCol the col of the center of rotation.
		 * @param rotation  the rotation in degrees.
		 * @return the minutia rotated around the given center.
		 */
		public static int[] applyRotation(int[] minutia, int centerRow, int centerCol, int rotation) {

			// calcul des formules nécessaires pour obtenir newRow, newCol et newOrientation
			int x = minutia[1] - centerCol;
			int y = centerRow - minutia[0];
			double newX = x * Math.cos(rotation*(Math.PI/180)) - y * Math.sin(rotation*(Math.PI/180));
			double newY = x * Math.sin(rotation*(Math.PI/180)) + y * Math.cos(rotation*(Math.PI/180));
			double newRow = centerRow - newY;
			double newCol = newX + centerCol;

			int Row = (int) Math.round(newRow);
			int Col = (int) Math.round(newCol);

			double newOrientation = ((minutia[2] + rotation) % 360);
			// transtypé en int : directement int ou rounded ????

			int[] paramètres = {Row, Col,(int) newOrientation};

			// retourne les coordonnées & l'orientation de la minutie après rotation
			return paramètres;
		}

		/**
		 * Applies the specified translation to the minutia.
		 *
		 * @param minutia        the original minutia.
		 * @param rowTranslation the translation along the rows.
		 * @param colTranslation the translation along the columns.
		 * @return the translated minutia.
		 */
		public static int[] applyTranslation ( int[] minutia, int rowTranslation, int colTranslation){
			// nouvelles coordonnées de la minutie après avoir apppliqué la translation
			int newRow = minutia[0] - rowTranslation;
			int newCol = minutia[1] - colTranslation;
			int newOrientation = minutia[2];

			int[] paramètres = {newRow, newCol, newOrientation};

			return paramètres;
		}

		/**
		 * Computes the row, column, and angle after applying a transformation
		 * (translation and rotation).
		 *
		 * @param minutia        the original minutia.
		 * @param centerCol      the column around which the point is rotated.
		 * @param centerRow      the row around which the point is rotated.
		 * @param rowTranslation the vertical translation.
		 * @param colTranslation the horizontal translation.
		 * @param rotation       the rotation.
		 * @return the transformed minutia.
		 */
		public static int[] applyTransformation ( int[] minutia, int centerRow, int centerCol, int rowTranslation,
		int colTranslation, int rotation){

			// crée un nouveau tableau pour Rotation (w/ tableau initial)
			int[] minutia2 = applyRotation(minutia, centerRow, centerCol, rotation);

			// crée un nouveau tableau pour Translation (w/ tableau de Rotation)
			int[] minutia3 = applyTranslation (minutia2, rowTranslation, colTranslation);

			return minutia3;
		}

		/**
		 * Computes the row, column, and angle after applying a transformation
		 * (translation and rotation) for each minutia in the given list.
		 *
		 * @param minutiae       the list of minutiae.
		 * @param centerCol      the column around which the point is rotated.
		 * @param centerRow      the row around which the point is rotated.
		 * @param rowTranslation the vertical translation.
		 * @param colTranslation the horizontal translation.
		 * @param rotation       the rotation.
		 * @return the list of transformed minutiae.
		 */
		public static List<int[]> applyTransformation (List < int[]>minutiae,int centerRow, int centerCol,
		int rowTranslation, int colTranslation, int rotation){

			//applique la transformation à une liste entière
			ArrayList <int[]> ListeMinutia = new ArrayList<int[]>();

			for(int[] i : minutiae) {
				ListeMinutia.add(applyTransformation(i, centerRow, centerCol, rowTranslation,
				colTranslation, rotation));
			}

			return ListeMinutia;
		}
		/**
		 * Counts the number of overlapping minutiae.
		 *
		 * @param minutiae1      the first set of minutiae.
		 * @param minutiae2      the second set of minutiae.
		 * @param maxDistance    the maximum distance between two minutiae to consider
		 *                       them as overlapping.
		 * @param maxOrientation the maximum difference of orientation between two
		 *                       minutiae to consider them as overlapping.
		 * @return the number of overlapping minutiae.
		 */
		public static int matchingMinutiaeCount (List < int[]>minutiae1, List < int[]>minutiae2,int maxDistance,
		int maxOrientation){

			// initialisation du nombre de match
			int matchingMinutiaes = 0;



			// deux boucles for qui vont comparer la distance et l'orientation entre les 2 minuties
			for(int i = 0; i < minutiae1.size(); ++i) {
				for(int j = 0; j < minutiae2.size(); ++j) {
					double distance = Math.sqrt((minutiae1.get(i)[0] - minutiae2.get(j)[0])*(minutiae1.get(i)[0]-minutiae2.get(j)[0])
							       + (minutiae1.get(i)[1] - minutiae2.get(j)[1])*(minutiae1.get(i)[1] - minutiae2.get(j)[1]));
					double orientationDiff = Math.abs(minutiae1.get(i)[2] - minutiae2.get(j)[2]);

					// si les condiitons sont respectées, alors les deux minuties match (+1)
					if(distance <= maxDistance && orientationDiff <= maxOrientation) {
						++matchingMinutiaes;
						break;
					}

				}

			}
			return matchingMinutiaes;
		}

		/**
		 * Compares the minutiae from two fingerprints.
		 *
		 * @param minutiae1 the list of minutiae of the first fingerprint.
		 * @param minutiae2 the list of minutiae of the second fingerprint.
		 * @return Returns <code>true</code> if they match and <code>false</code>
		 *         otherwise.
		 */
		public static boolean match (List < int[]>minutiae1, List < int[]>minutiae2){

			// initialisation des variables de base
			int maxDistance = DISTANCE_THRESHOLD;
			int maxOrientation = ORIENTATION_DISTANCE;

			// deux boucles for qui vont prendre chaque élément de m1, auxquels on applique la Transformation
			// et que l'on va comparer à chaque minutie de m2
			for (int i = 0; i < minutiae2.size(); ++i) {
				for (int j = 0; j < minutiae1.size(); ++j) {

					//déclaration des variables utilisées
					// minutiae1[i] : row élément 0 et col élément 1
					int centerRow = minutiae2.get(i)[0];
					int centerCol = minutiae2.get(i)[1];
					int rowTranslation = minutiae2.get(i)[0] - minutiae1.get(j)[0];
					int colTranslation = minutiae2.get(i)[1] - minutiae1.get(j)[1];
					int rotation = minutiae2.get(i)[2] - minutiae1.get(j)[2];

					// boucle for qui applique la transformation de chaque minutie de la liste 1
					for(int g = rotation - MATCH_ANGLE_OFFSET; g <= rotation + MATCH_ANGLE_OFFSET; ++g ) {
						int matches = matchingMinutiaeCount(applyTransformation(minutiae2, centerRow, centerCol, rowTranslation, colTranslation, g),
											minutiae1, maxDistance, maxOrientation);
						boolean same = (matches >= FOUND_THRESHOLD);

						if(same) {
							System.out.print("Match count "+matches+" ");
							return true;
						}
					}
				}
			}

			// 2 boucles for qui testent chaque minutie
			// on fait apply transformation pour chaque minutie de la liste 1, vers chaque minutie de la liste 2
			// on compte le nombre de matching des autres minuties dans chaque cas
			// if nombre trouvé avec matching == distance_thresold, alors true

			return false;
		}

}
