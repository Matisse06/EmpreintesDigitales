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
		if ((row < 0 || row >= image.length) || (col < 0 || col>= image[0].length)) {
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
			if (neighbours[i] == true) {
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

			if (neighbours[i] == false && neighbours[j] == true) {
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


				if (pixelNoir) {
					System.out.println("case " + i + ", " + j);
					System.out.println("blackN : " + blackNeighbours(neighbours));
					System.out.println("nb transitions : " + transitions(getNeighbours(imageTest, i, j))); }


				if ((pixelNoir) && (neighbours != null) &&
						(blackNeighbours(neighbours)) >= 2 && (blackNeighbours(neighbours) <= 6) &&
						(transitions(neighbours) == 1)) {

					// test step 1

					if (step == 0) {

						if (((getNeighbours(imageTest, i, j)[0] == false) ||(getNeighbours(imageTest, i, j)[2] == false) || (getNeighbours(imageTest, i, j)[4] == false)) &&
								((getNeighbours(imageTest, i, j)[2] == false) ||(getNeighbours(imageTest, i, j)[4] == false) || (getNeighbours(imageTest, i, j)[6] == false))) {
							imageTest4[i][j] = false;
							System.out.println("deleted !");
						}


						// test step 2

					} else if (step == 1) {


						if (((getNeighbours(imageTest, i, j)[0] == false) ||(getNeighbours(imageTest, i, j)[2] == false) || (getNeighbours(imageTest, i, j)[6] == false)) &&
								((getNeighbours(imageTest, i, j)[0] == false) ||(getNeighbours(imageTest, i, j)[4] == false) || (getNeighbours(imageTest, i, j)[6] == false))) {
							imageTest4[i][j] = false;
							System.out.println("deleted !");
						}


					}

					System.out.println();
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
 /* public static boolean[][] thin(boolean[][] image) {
	  //TODO implement
	  return null;
  }

  /**
   * Computes all pixels that are connected to the pixel at coordinate
   * <code>(row, col)</code> and within the given distance of the pixel.
   *
   * @param image    array containing each pixel's boolean value.
   * @param row      the first coordinate of the pixel of interest.
   * @param col      the second coordinate of the pixel of interest.
   * @param distance the maximum distance at which a pixel is considered.
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
	public int somme(int a, int b){
		return 0;
	}

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
				System.out.println("Pix ligne" + i);
				for (int j = 0; j < image[0].length; ++j) {
					System.out.println("Pix col" + j);
					if (image[i][j] == false) {
						//pixel blanc, go pixel suivant
						System.out.println("pixel blanc");
					} else {
						System.out.println("  ");
						System.out.println("  ");
						Main.printArray(connected);
						System.out.println("  ");
						if ((i < row - distance || i > row + distance) || (j < col - distance || j > col + distance)) {
							// trop loin de la minutiae, go suivant
							System.out.println("pixel out of bounds");
						} else {
							if (connected[i][j] == true) {
								// si pixel noir deja dans tableau connected, go pixel suivant
								System.out.println("pixel noir deja dans tab");
							} else {


								boolean[] neighbours = getNeighbours(image, i, j);
								for (int g = 0; g < neighbours.length; ++g) {
									if (neighbours[g] == true) {
										int rowPixNoir = rowNeighbour(g, i, j);                // determiner position pixel noir voisin
										int colPixNoir = colNeighbour(g, i, j);
										//System.out.println("coordonnees pixel noir connecte" +rowPixNoir+ " " +colPixNoir);
										if (connected[rowPixNoir][colPixNoir] == true) {        // si le pixel voisin est lie a minutiae good!
											connected[i][j] = true;
											System.out.println("add lie dans tab");
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

			if (somme_x2 == 0) {
				slope = Double.POSITIVE_INFINITY;
				return slope;
			} else {

				if (somme_x2 >= somme_y2) {
					slope = somme_xy / somme_x2;
				} else {
					slope = somme_y2 / somme_xy;
				}

				return slope;
			}
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

			// calcul de l'angle formé par le vecteur direction
			double angle = Math.atan(slope);


			// variable qui compte le nombre de pixels au dessus et en dessous
			int pixelUp = 0;
			int pixelDown = 0;

			// boucle qui verifie que les conditions soient respectées

			for (int row = 0; row < connectedPixels.length; ++row) {
				for (int col = 0; col < connectedPixels[row].length; ++col) {

					// x et y exprimés en fonction de l'origine du repère, centré sur la minutie
					int x = col - colM;
					int y = rowM - row;

					boolean pixelNoir = connectedPixels[row][col];

					if (pixelNoir && (y >= ((-1 / slope) * x))) {
						pixelUp += 1;
					} else if (pixelNoir && (y < ((-1 / slope) * x))) {
						pixelDown += 1;
					}
				}
			}

			// conditions qui déterminent l'ajout potentiel de PI à l'angle calculé au début
			if ((slope == Double.POSITIVE_INFINITY && (pixelUp > pixelDown))) {
				// on est d'accord, le "au dessus de la minutie" c'est bien de la perpendiculaire ???
				angle = Math.PI / 2;
			} else if ((slope == Double.POSITIVE_INFINITY && (pixelUp < pixelDown))) {
				angle = -Math.PI / 2;
			} else if ((angle > 0 && (pixelDown > pixelUp)) || (angle < 0 && (pixelDown < pixelUp))) {
				// vérifier si c'est strictement sup / inf ou pas !!!!! si non = what if == ????
				angle += Math.PI;
			}

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
			//
			// angle en radian
			//double angleR = computeAngle(computeSlope(connectedPixels(image, row, col, distance)));
			//angleR = Math.toDegrees(angleR);

			// angle en dégré (rounded)
			//int angleD = Math.round(angleD);
			// initialisé à double mais comment le rendre int ?

			// if (angleD < 0) {
			//	  angleD += 360;
			//}

			//return angleD;
			// comment savoir les coordonnées de la minutie ??

			//
			return 0;
		}

		/**
		 * Extracts the minutiae from a thinned image.
		 *
		 * @param image array containing each pixel's boolean value.
		 * @return The list of all minutiae. A minutia is represented by an array where
		 *         the first element is the row, the second is column, and the third is
		 *         the angle in degrees.
		 * @see #thin(boolean[][])
		 */
		public static List<int[]> extract ( boolean[][] image){
			//TODO implement
			return null;
		}

		/**
		 *
		 *
		 *
		 *
		 *
		 *
		 *
		 * Applies the specified rotation to the minutia.
		 *
		 * @param minutia   the original minutia.
		 * @param centerRow the row of the center of rotation.
		 * @param centerCol the col of the center of rotation.
		 * @param rotation  the rotation in degrees.
		 * @return the minutia rotated around the given center.
		 */
		public static int[] applyRotation ( int[] minutia, int centerRow, int centerCol, int rotation){
			//TODO implement
			return null;
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
			//TODO implement
			return null;
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
			//TODO implement
			return null;
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
		int rowTranslation,
		int colTranslation, int rotation){
			//TODO implement
			return null;
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
			//TODO implement
			return 0;
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
			//TODO implement
			return false;
		}

}
