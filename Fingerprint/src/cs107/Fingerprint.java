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
    public static boolean[] getNeighbours(boolean[][] image, int row, int col) {
        assert (image != null);
        assert (row >= 0 && row < image.length && col >= 0 && col < image[0].length);

        boolean[] neighbours = new boolean[8];
        neighbours[0] = color(image, row - 1, col);
        for (int i = -1; i < 2; ++i) {
            neighbours[i + 2] = color(image, row + i, col + 1);
        }
        neighbours[4] = color(image, row + 1, col);
        neighbours[5] = color(image, row + 1, col - 1);
        neighbours[6] = color(image, row, col - 1);
        neighbours[7] = color(image, row - 1, col - 1);

        return neighbours;
    }

    /**
     * Retunrs the boolean value of the pixel (false if outside of image). Method used in getNeighbours.
     * @param image array containing each pixel's boolean value.
     * @param row of pixel.
     * @param col of pixel.
     * @return boolean value of pixel.
     */
    public static boolean color(boolean[][] image, int row, int col) {
        assert (image != null);

        boolean couleur = false;
        if ((row >= image.length || row < 0) || (col >= image[0].length || col < 0)) {
            couleur = false;
        } else {
            couleur = image[row][col];
        }
        return couleur;
    }

    /**
     * Computes the number of black (<code>true</code>) pixels among the neighbours
     * of a pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   the convention described in
     *                   {@link #getNeighbours(boolean[][], int, int)}.
     * @variable black   number of black neighbours
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
     * <p>
     * pixel.
     *
     * @param neighbours array containing each pixel value. The array must respect
     *                   <p>
     *                   the convention described in
     *                   <p>
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
                ++numberOfTransitions;
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
     * <p>
     * otherwise.
     */

    public static boolean identical(boolean[][] image1, boolean[][] image2) {
        assert ((image1 != null) && (image2 != null));

        if ((image1.length != image2.length) || (image1[0].length != image2[0].length))
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
     * @param step    the step to apply, Step 0 or Step 1.
     * @variable imageTest4   copy of imageTest that endures the modifications
     * @return A new array containing each pixel's value after the step.
     */

    public static boolean[][] thinningStep(boolean[][] imageTest, int step) {
        assert ((step == 0 || step == 1) && (imageTest != null));

        boolean[][] imageTest4 = new boolean[imageTest.length][imageTest[0].length];

        for (int i = 0; i < imageTest.length; ++i) {
            for (int j = 0; j < imageTest[i].length; ++j) {
                imageTest4[i][j] = imageTest[i][j];
            }
        }

        for (int i = 0; i < imageTest.length; ++i) {
            for (int j = 0; j < imageTest[i].length; ++j) {

                boolean pixelNoir = imageTest[i][j];
                boolean[] neighbours = getNeighbours(imageTest, i, j);

                if ((pixelNoir) && (neighbours != null) &&
                        (blackNeighbours(neighbours)) >= 2 && (blackNeighbours(neighbours) <= 6) &&
                        (transitions(neighbours) == 1)) {

                    if (step == 0) {

                        if (((getNeighbours(imageTest, i, j)[0] == false) || (getNeighbours(imageTest, i, j)[2] == false) || (getNeighbours(imageTest, i, j)[4] == false)) &&
                                ((getNeighbours(imageTest, i, j)[2] == false) || (getNeighbours(imageTest, i, j)[4] == false) || (getNeighbours(imageTest, i, j)[6] == false))) {
                            imageTest4[i][j] = false;
                        }

                    } else if (step == 1) {

                        if (((getNeighbours(imageTest, i, j)[0] == false) || (getNeighbours(imageTest, i, j)[2] == false) || (getNeighbours(imageTest, i, j)[6] == false)) &&
                                ((getNeighbours(imageTest, i, j)[0] == false) || (getNeighbours(imageTest, i, j)[4] == false) || (getNeighbours(imageTest, i, j)[6] == false))) {
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
     *
     * @param image array containing each pixel's boolean value.
     * @variable imageTest1    copy of image
     * @variable imageTest2    tab on which step 0 is applied (w/ imageTest1 as input)
     * @variable imageTest3    tab on which step 1 is applied (w/ imageTest2 as input)
     * @return array containing the boolean value of each pixel of the image after
     * applying the thinning algorithm.
     */

    public static boolean[][] thin(boolean[][] image) {
        assert (image != null);

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
     * Returns the row (coordinate in image) of pixel depending on its position in array returned by getNeighbours.
     * @param i position of the pixel in array returned by getNeighbours.
     * @param row of pixel on which getNeighbours was used.
     * @param col of pixel on which getNeighbours was used.
     * @return int value corresponding to the row of pixel in image.
     */
    public static int rowNeighbour(int i, int row, int col) {

        if (i == 0 || i == 1 || i == 7) {
            return row - 1;
        }
        if (i == 2 || i == 6) {
            return row;
        }
        if (i == 3 || i == 4 || i == 5) {
            return row + 1;
        }
        return -1;
    }

    /**
     * Returns the col (coordinate in image) of pixel depending on its position in array returned by getNeighbours.
     * @param i position of the pixel in array returned by getNeighbours.
     * @param row of pixel on which getNeighbours was used.
     * @param col of pixel on which getNeighbours was used.
     * @return int value corresponding to the col of pixel in image.
     */
    public static int colNeighbour(int i, int row, int col) { // determine la colonne du pixel voisin de getNeighbours
        if (i == 0 || i == 4) {
            return col;
        }
        if (i == 1 || i == 2 || i == 3) {
            return col + 1;
        }
        if (i == 5 || i == 6 || i == 7) {
            return col - 1;
        }
        return -1;
    }

   /**
     * Computes all pixels that are connected to the pixel at coordinate
     * <code>(row, col)</code> and within the given distance of the pixel.
     *
     * @param //image    array containing each pixel's boolean value.
     * @param row        the first coordinate of the pixel of interest.
     * @param col        the second coordinate of the pixel of interest.
     * @param //distance the maximum distance at which a pixel is considered.
     * @return An array where <code>true</code> means that the pixel is within
     * <code>distance</code> and connected to the pixel at
     * <code>(row, col)</code>.
     */
    public static boolean[][] connectedPixels(boolean[][] image, int row, int col, int distance) {
       	assert (image != null);
        assert (row >= 0 && row < image.length && col >=0 && col < image[0].length);

        boolean[][] connected = new boolean[image.length][image[0].length];
        connected[row][col] = true;

        boolean finit = false;
        while (finit == false) {
            finit = true;
            for (int i = 0; i < image.length; i++) {
                for (int j = 0; j < image[0].length; ++j) {
                    if (image[i][j] == false) {
                        //if pixel is white, go to the next pixel
                    } else {
                        if ((i < row - distance || i > row + distance) || (j < col - distance || j > col + distance)) {
                            // if pixel too far from the minutiae, go next pixel
                        } else {
                            if (connected[i][j] == true) {
                                // if pixel already in array connected, go to the following pixel
                            } else {
                                boolean[] neighbours = getNeighbours(image, i, j);
                                for (int g = 0; g < neighbours.length; ++g) {
                                    if (neighbours[g] == true) {
                                        int rowPixNoir = rowNeighbour(g, i, j);
                                        int colPixNoir = colNeighbour(g, i, j);
                                        if (connected[rowPixNoir][colPixNoir] == true) {
                                            connected[i][j] = true;
                                            finit = false;
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
     * @param //row           the row of the minutia.
     * @param //col           the col of the minutia.
     * @return the slope.
     */
    public static double computeSlope(boolean[][] connectedPixels, int rowM, int colM) {
        assert (connectedPixels != null);
        assert(rowM >= 0 && rowM < connectedPixels.length && colM >= 0 && colM < connectedPixels[0].length);

        double slope = 0.;

        double somme_x2 = 0;
        double somme_y2 = 0;
        double somme_xy = 0;


        for (int row = 0; row < connectedPixels.length; ++row) {
            for (int col = 0; col < connectedPixels[row].length; ++col) {

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
     * @param //row           the row of the minutia.
     * @param //col           the col of the minutia.
     * @param slope           the slope as returned by
     *                        {@link #computeSlope(boolean[][], int, int)}.
     * @variable angle        angle (rad) initialement calculé avec atan(slope)
     * @variable pixelUp      nombre de pixels situés au dessus de la slope
     * @variable pixelDown    nombre de pixels situés en dessous de la slope
     * @return the orientation of the minutia in radians.
     */


    public static double computeAngle(boolean[][] connectedPixels, int rowM, int colM, double slope) {
        assert (connectedPixels != null);
        assert(rowM >= 0 && rowM < connectedPixels.length && colM >= 0 && colM < connectedPixels[0].length);


        double angle = Math.atan(slope);
        int pixelUp = 0;
        int pixelDown = 0;

        for (int row = 0; row < connectedPixels.length; ++row) {
            for (int col = 0; col < connectedPixels[row].length; ++col) {

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

        if ((slope == Double.POSITIVE_INFINITY && (pixelUp > pixelDown))) {
            angle = Math.PI / 2;

        } else if ((slope == Double.POSITIVE_INFINITY && (pixelUp < pixelDown))) {
            angle = -Math.PI / 2;

        } else if ((angle > 0 && (pixelDown > pixelUp)) || (angle < 0 && (pixelDown <= pixelUp))) {
            angle += Math.PI;
        } else if ((angle == 0 && (pixelDown > pixelUp))) {
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
     * @variable angleR  computed angle in radian, transformed into degrees
     * @variable angleD  computed angle in degrees, rounded
     *
     * @return The orientation in degrees.
     */
    public static int computeOrientation(boolean[][] image, int row, int col, int distance) {
        assert (image != null);
        assert (row >= 0 && row < image.length && col >= 0 && col < image[0].length);

        boolean[][] neighbours = connectedPixels(image, row, col, distance);
        double slope = computeSlope(neighbours, row, col);
        double angleR = computeAngle(neighbours, row, col, slope);

        angleR = Math.toDegrees(angleR);

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
     * the first element is the row, the second is column, and the third is
     * the angle in degrees.
     * @see #thin(boolean[][])
     */
    public static List<int[]> extract(boolean[][] image) {
        assert (image != null);

        ArrayList<int[]> listeMinutie = new ArrayList<int[]>();

        for (int i = 1; i < image.length - 1; ++i) {
            for (int j = 1; j < image[0].length - 1; ++j) {
                if (image[i][j]) {
                    int transi = transitions(getNeighbours(image, i, j));
                    if (transi == 1 || transi == 3) {
                        int[] minutie = new int[3];
                        minutie[0] = i;
                        minutie[1] = j;
                        minutie[2] = computeOrientation(image, i, j, ORIENTATION_DISTANCE);
                        listeMinutie.add(minutie);
                    }
                }
            }
        }
        return listeMinutie;    // returns list of minuatiaes( their location + orientation)
    }

    /**
     * Applies the specified rotation to the minutia.
     *
     * @param minutia   the original minutia.
     * @param centerRow the row of the center of rotation.
     * @param centerCol the col of the center of rotation.
     * @param rotation  the rotation in degrees.
     * @variable paramètres    parameters of the rotated minutia
     * @return the minutia rotated around the given center.
     */
    public static int[] applyRotation(int[] minutia, int centerRow, int centerCol, int rotation) {
        assert (minutia != null);
        assert (centerRow >= 0 && centerCol >= 0);

        int x = minutia[1] - centerCol;
        int y = centerRow - minutia[0];
        double newX = x * Math.cos(rotation * (Math.PI / 180)) - y * Math.sin(rotation * (Math.PI / 180));
        double newY = x * Math.sin(rotation * (Math.PI / 180)) + y * Math.cos(rotation * (Math.PI / 180));
        double newRow = centerRow - newY;
        double newCol = newX + centerCol;

        int Row = (int) Math.round(newRow);
        int Col = (int) Math.round(newCol);

        double newOrientation = ((minutia[2] + rotation) % 360);

        int[] paramètres = {Row, Col, (int) newOrientation};

        return paramètres;
    }

    /**
     * Applies the specified translation to the minutia.
     *
     * @param minutia        the original minutia.
     * @param rowTranslation the translation along the rows.
     * @param colTranslation the translation along the columns.
     * @variable paramètres1 parameters of translated minutia
     * @return the translated minutia.
     */
    public static int[] applyTranslation(int[] minutia, int rowTranslation, int colTranslation) {
        assert (minutia != null);
        assert (rowTranslation >= 0 && colTranslation >= 0);

        int newRow = minutia[0] - rowTranslation;
        int newCol = minutia[1] - colTranslation;
        int newOrientation = minutia[2];

        int[] paramètres1 = {newRow, newCol, newOrientation};

        return paramètres1;
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
     * @variable minutiaR    minutia after Rotation
     * @variable minutiaT    minutia after Translation & Rotation
     * @return the transformed minutia.
     */
    public static int[] applyTransformation(int[] minutia, int centerRow, int centerCol, int rowTranslation,
                                            int colTranslation, int rotation) {

        assert (minutia != null);
        assert (centerRow >= 0 && centerCol >=0 && rowTranslation >= 0 && colTranslation >= 0);

        int[] minutiaR = applyRotation(minutia, centerRow, centerCol, rotation);

        int[] minutiaT = applyTranslation(minutiaR, rowTranslation, colTranslation);

        return minutiaT;
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
    public static List<int[]> applyTransformation(List<int[]> minutiae, int centerRow, int centerCol,
                                                  int rowTranslation, int colTranslation, int rotation) {

        assert (minutiae != null);
        assert (centerRow >= 0 && centerCol >=0 && rowTranslation >= 0 && colTranslation >= 0);

        ArrayList<int[]> ListeMinutia = new ArrayList<int[]>();

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
    public static int matchingMinutiaeCount(List<int[]> minutiae1, List<int[]> minutiae2, int maxDistance,
                                            int maxOrientation) {

       assert (minutiae1 != null && minutiae2 != null);

        int matchingMinutiaes = 0;
        boolean[] matched = new boolean[minutiae2.size()];

        for (int i = 0; i < minutiae1.size(); ++i) {
            for (int j = 0; j < minutiae2.size(); ++j) {
                double distance = Math.sqrt((minutiae1.get(i)[0] - minutiae2.get(j)[0]) * (minutiae1.get(i)[0] - minutiae2.get(j)[0])
                        + (minutiae1.get(i)[1] - minutiae2.get(j)[1]) * (minutiae1.get(i)[1] - minutiae2.get(j)[1]));
                double orientationDiff = Math.abs(minutiae1.get(i)[2] - minutiae2.get(j)[2]);

                if (distance <= maxDistance && orientationDiff <= maxOrientation && !matched[j]) {
                    ++matchingMinutiaes;
                    matched[j] = true;
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
     * @variable matches    number of match between the two fingerprints
     * @return Returns <code>true</code> if they match and <code>false</code>
     * otherwise.
     */
    public static boolean match(List<int[]> minutiae1, List<int[]> minutiae2) {
        assert (minutiae1 != null && minutiae2 != null);

        int maxDistance = DISTANCE_THRESHOLD;
        int maxOrientation = ORIENTATION_THRESHOLD;

        for (int i = 0; i < minutiae2.size(); ++i) {
            for (int j = 0; j < minutiae1.size(); ++j) {

                int centerRow = minutiae2.get(i)[0];
                int centerCol = minutiae2.get(i)[1];
                int rowTranslation = minutiae2.get(i)[0] - minutiae1.get(j)[0];
                int colTranslation = minutiae2.get(i)[1] - minutiae1.get(j)[1];
                int rotation = minutiae2.get(i)[2] - minutiae1.get(j)[2];


                for (int g = rotation - MATCH_ANGLE_OFFSET; g <= rotation + MATCH_ANGLE_OFFSET; ++g) {
                    int matches = matchingMinutiaeCount(applyTransformation(minutiae2, centerRow, centerCol, rowTranslation, colTranslation, g),
                            minutiae1, maxDistance, maxOrientation);

                    boolean same = (matches >= FOUND_THRESHOLD);
                    if (same) {
                        return true;
                    }
                }

            }
        }
        return false;
    }

}
