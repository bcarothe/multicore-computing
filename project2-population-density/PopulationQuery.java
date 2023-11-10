
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.RecursiveTask;
import java.util.concurrent.RecursiveAction;

// ==================================================
// VERSION 2
// ==================================================

class V2_GetCorners extends RecursiveAction {
    // static final int SEQUENTIAL_THRESHOLD = 5000;

    // Member variables
    int low;
    int high;
    CensusData census_data;
    float[] corners;

    // Constructor
    V2_GetCorners(CensusData census_d, int lo, int hi) {
        census_data = census_d;
        low = lo;
        high = hi;
    }

    // Override of the RecursiveAction compute() function
    protected void compute() {
        corners = new float[4];

        if (high - low == 1) {
            // float sum = 0;
            // for (int i = low; i < high; ++i)
            // sum += array[i];
            // return sum;

            float latitude = census_data.data[low].latitude;
            float longitude = census_data.data[low].longitude;
            corners = new float[] { longitude, latitude, longitude, latitude };
        } else {
            // int mid = low + (high - low) / 2;
            // V2_GetCorners left = new V2_GetCorners(array, low, mid);
            // V2_GetCorners right = new V2_GetCorners(array, mid, high);
            // left.fork();
            // float rightAns = right.compute();
            // float leftAns = left.join();
            // return leftAns + rightAns;

            V2_GetCorners left = new V2_GetCorners(census_data, low, (low + high) / 2);
            V2_GetCorners right = new V2_GetCorners(census_data, (low + high) / 2, high);

            left.fork();
            right.compute();
            left.join();

            // Iterate over each corner values
            for (int i = 0; i < corners.length; i++) {
                // First 2 corner values
                if (i < corners.length / 2)
                    if (left.corners[i] < right.corners[i]) {
                        corners[i] = left.corners[i];
                    } else {
                        corners[i] = right.corners[i];
                    }
                // Last 2 corner values
                else {
                    if (left.corners[i] > right.corners[i]) {
                        corners[i] = left.corners[i];
                    } else {
                        corners[i] = right.corners[i];
                    }
                }
            }
        }
    }

    // Invokes the fork-join process and returns the calculated corners
    static float[] getCorners(CensusData census_data) {
        V2_GetCorners v2_getCorners = new V2_GetCorners(census_data, 0, census_data.data_size);
        ForkJoinPool.commonPool().invoke(v2_getCorners);
        return v2_getCorners.corners;
    }

}

class V2_GetQueryResults extends RecursiveAction {
    // Member variables
    float[] results = new float[2];
    float subset_population = 0;
    float total_population = 0;
    int low;
    int high;
    CensusData census_data;
    float[] corners;
    int row_count, column_count;
    int west;
    int south;
    int east;
    int north;

    // Constructor
    V2_GetQueryResults(CensusData census_d, int rc, int cc, float[] corns, int w, int s, int e, int n, int lo, int hi) {
        census_data = census_d;
        row_count = rc;
        column_count = cc;
        corners = corns;
        west = w;
        south = s;
        east = e;
        north = n;
        low = lo;
        high = hi;
    }

    // Override of the RecursiveAction compute() function
    protected void compute() {
        if (high - low == 1) {
            // Add the total population
            total_population += census_data.data[low].population;

            // Get block coords
            int blockX = (int) ((census_data.data[low].longitude - corners[0]) / ((corners[2] - corners[0]) / row_count)
                    + 1);
            int blockY = (int) ((census_data.data[low].latitude - corners[1])
                    / ((corners[3] - corners[1]) / column_count) + 1);

            // Add population to subset population if the block is in the query rectangle
            if (blockX >= west && blockX < east + 1 && blockY >= south && blockY < north + 1) {
                subset_population += census_data.data[low].population;
            }
        } else {
            V2_GetQueryResults left = new V2_GetQueryResults(census_data, row_count, column_count, corners, west, south,
                    east, north, low, (low + high) / 2);
            V2_GetQueryResults right = new V2_GetQueryResults(census_data, row_count, column_count, corners, west,
                    south, east, north, (low + high) / 2, high);

            left.fork();
            right.compute();
            left.join();

            // Add left total and subset pupulation results
            subset_population = left.subset_population + right.subset_population;
            total_population = left.total_population + right.total_population;
        }
    }

    // Invokes the fork-join process and returns the calculated corners
    static float[] getQueryResults(CensusData census_data, int row_count, int column_count, float[] corners, int west,
            int south, int east, int north) {
        V2_GetQueryResults V2_getQueryResults = new V2_GetQueryResults(census_data, row_count, column_count, corners,
                west, south, east, north, 0, census_data.data_size);
        ForkJoinPool.commonPool().invoke(V2_getQueryResults);
        V2_getQueryResults.results[0] = V2_getQueryResults.subset_population;
        V2_getQueryResults.results[1] = V2_getQueryResults.subset_population / V2_getQueryResults.total_population;
        return V2_getQueryResults.results;
    }

}

// //////////////////////////////////////////////////

// ==================================================
// VERSION 4
// ==================================================

class V4_CreateGrid extends RecursiveAction {
    // Member variables
    static final int SEQUENTIAL_THRESHOLD = 5000;
    float subset_population = 0;
    float total_population = 0;
    int low;
    int high;
    CensusData census_data;
    float[] corners;
    int row_count, column_count;
    int west;
    int south;
    int east;
    int north;
    int[][] grid;

    // Constructor
    public V4_CreateGrid(CensusData census_d, int rc, int cc, float[] corns, int w, int s, int e, int n, int lo,
            int hi) {
        census_data = census_d;
        row_count = rc;
        column_count = cc;
        corners = corns;
        west = w;
        south = s;
        east = e;
        north = n;
        low = lo;
        high = hi;
        grid = new int[row_count][column_count];
    }

    // Override of the RecursiveAction compute() function
    public void compute() {
        if (high - low < SEQUENTIAL_THRESHOLD) {
            for (int i = low; i < high; i++) {
                // Get block coords from each CensusData point
                int blockX = (int) ((census_data.data[i].longitude - corners[0])
                        / ((corners[2] - corners[0]) / row_count));
                int blockY = (int) ((census_data.data[i].latitude - corners[1])
                        / ((corners[3] - corners[1]) / column_count));

                // Convert block coords to array indices
                int blockX_i = blockX - 1;
                int blockY_i = blockY - 1;

                // Use array indices at edges to avoid out of bounds
                if (blockX == row_count && blockY == column_count) { // Top right corner
                    grid[blockX_i][blockY_i] += census_data.data[i].population;
                }
                // East edge
                else if (blockX == row_count) {
                    grid[blockX_i][blockY] += census_data.data[i].population;
                }
                // north edge
                else if (blockY == column_count) {
                    grid[blockX][blockY_i] += census_data.data[i].population;
                }
                // Not on east or north edge
                else {
                    grid[blockX][blockY] += census_data.data[i].population;
                }
            }
        } else {
            V4_CreateGrid left = new V4_CreateGrid(census_data, row_count, column_count, corners, west, south, east,
                    north, low, (low + high) / 2);
            V4_CreateGrid right = new V4_CreateGrid(census_data, row_count, column_count, corners, west, south, east,
                    north, (low + high) / 2, high);

            left.fork();
            right.compute();
            left.join();

            grid = V4_CombineGrid.combineGrid(left.grid, right.grid);
        }
    }
}

class V4_CombineGrid extends RecursiveAction {
    public static final int SEQUENTIAL_THRESHOLD = 16;
    int[][] grid_left;
    int[][] grid_right;
    int[][] grid_combined;
    int row_low;
    int row_high;
    int column_low;
    int column_high;

    public V4_CombineGrid(int[][] grid_l, int[][] grid_r, int[][] grid_c, int row_l, int row_h, int col_l, int col_h) {
        grid_left = grid_l;
        grid_right = grid_r;
        grid_combined = grid_c;
        row_low = row_l;
        row_high = row_h;
        column_low = col_l;
        column_high = col_h;
    }

    public void compute() {
        if (row_high - row_low < SEQUENTIAL_THRESHOLD || column_high - column_low < SEQUENTIAL_THRESHOLD) {
            for (int i = row_low; i < row_high; i++) {
                for (int j = column_low; j < column_high; j++) {
                    grid_combined[i][j] += grid_left[i][j] + grid_right[i][j];
                }
            }
        } else {
            V4_CombineGrid bottom_left = new V4_CombineGrid(grid_left, grid_right, grid_combined, row_low,
                    (row_low + row_high) / 2, column_low, (column_low + column_high) / 2);
            V4_CombineGrid bottom_right = new V4_CombineGrid(grid_left, grid_right, grid_combined, row_low,
                    (row_low + row_high) / 2, (column_low + column_high) / 2, column_high);
            V4_CombineGrid top_left = new V4_CombineGrid(grid_left, grid_right, grid_combined, (row_low + row_high) / 2,
                    row_high, column_low, (column_low + column_high) / 2);
            V4_CombineGrid top_right = new V4_CombineGrid(grid_left, grid_right, grid_combined,
                    (row_low + row_high) / 2, row_high, (column_low + column_high) / 2, column_high);

            bottom_left.fork();
            bottom_right.compute();
            top_left.compute();
            top_right.compute();
            bottom_left.join();
        }
    }

    public static int[][] combineGrid(int[][] grid_left, int[][] grid_right) {
        int sizeX = grid_left.length;
        int sizeY = grid_left[0].length;
        int[][] grid_combined = new int[sizeX][sizeY];

        V4_CombineGrid V4_combineGrid = new V4_CombineGrid(grid_left, grid_right, grid_combined, 0, sizeX, 0, sizeY);
        ForkJoinPool.commonPool().invoke(V4_combineGrid);

        return V4_combineGrid.grid_combined;
    }
}

class V4_GetQueryResults {
    // Member variables
    float[] results = new float[2];
    float subset_population = 0;
    float total_population = 0;
    int low;
    int high;
    CensusData census_data;
    float[] corners;
    int row_count, column_count;
    int west;
    int south;
    int east;
    int north;

    // Constructor
    V4_GetQueryResults(CensusData census_d, int rc, int cc, float[] corns, int w, int s, int e, int n, int lo, int hi) {
        census_data = census_d;
        row_count = rc;
        column_count = cc;
        corners = corns;
        west = w;
        south = s;
        east = e;
        north = n;
        low = lo;
        high = hi;
    }

    public static float[] getQueryResults(CensusData census_data, int row_count, int column_count, float[] corners,
            int west, int south, int east, int north) {
        //////////////////
        // Create grid
        V4_CreateGrid V4_createGrid = new V4_CreateGrid(census_data, row_count, column_count, corners, west, south,
                east, north, 0, census_data.data_size);
        ForkJoinPool.commonPool().invoke(V4_createGrid);
        int[][] grid = V4_createGrid.grid;

        ////////////////////////
        // Calculate results
        float[] result = new float[2];
        int bottom_right;
        int top_right;
        int bottom_left;
        int top_left;

        // Bottom-right of query rectangle
        bottom_right = grid[east - 1][north - 1];

        // Top-right of query rectangle
        if (south - 2 < 0) {
            top_right = 0;
        } else {
            top_right = grid[east - 1][south - 2];
        }

        // Bottom-left of query rectangle
        if (west - 2 < 0) {
            bottom_left = 0;
        } else {
            bottom_left = grid[west - 2][north - 1];
        }

        // Top-left of query rectangle
        if (south - 2 < 0 || west - 2 < 0) {
            top_left = 0;
        } else {
            top_left = grid[west - 2][south - 2];
        }

        // Store results
        result[0] = bottom_right - top_right - bottom_left + top_left; // Population in query rectangle
        result[1] = result[0] / (grid[row_count - 1][column_count - 1]); // Percentage of population in query rectangle

        return result;
    }
}

// //////////////////////////////////////////////////

// ==================================================
// VERSION 5
// ==================================================

class V5_GetQueryResults {
    // Member variables
    CensusData census_data;
    int row_count;
    int column_count;
    float[] corners;
    int west;
    int south;
    int east;
    int north;
    int thread_count;

    int[][] grid;
    int[][] locks;

    // Constructor
    public V5_GetQueryResults(CensusData census_d, int row_c, int column_c, float[] corns, int w, int s, int e, int n,
            int thread_c) {
        census_data = census_d;
        row_count = row_c;
        column_count = column_c;
        corners = corns;
        west = w;
        south = s;
        east = e;
        north = n;
        thread_count = thread_c;
        grid = new int[row_count][column_count];
        locks = new int[row_count][column_count];

        // fill the locks with all 0s
        for (int x = 0; x < locks.length; x++) {
            for (int y = 0; y < locks[0].length; y++)
                locks[x][y] = 0;
        }
    }

    public float[] V5_getQueryResults() {
        // Hold the provided number of threads
        V5_CreateGrid grids[] = new V5_CreateGrid[thread_count];

        // Split and start the threads
        for (int i = 0; i < thread_count - 1; i++) {
            int l = i * (census_data.data_size / thread_count);
            int h = (i + 1) * (census_data.data_size / thread_count);

            grids[i] = new V5_CreateGrid(census_data, row_count, column_count, corners, west, south, east, north, grid,
                    locks, l, h);
            grids[i].start();
        }

        int l = (thread_count - 1) * (census_data.data_size / thread_count);
        int h = census_data.data_size;

        grids[thread_count - 1] = new V5_CreateGrid(census_data, row_count, column_count, corners, west, south, east,
                north, grid, locks, l, h);
        grids[thread_count - 1].start();

        // Join the threads
        for (int i = 0; i < thread_count; i++) {
            try {
                grids[i].join();
            } catch (InterruptedException exception) {
                throw new RuntimeException(exception);
            }
        }

        // Process the grid
        grid = PopulationQuery.v3_ProcessGrid(grid);

        ///////////////////////
        // Calculate results
        float[] result = new float[2];
        int bottom_right;
        int top_right;
        int bottom_left;
        int top_left;

        // Bottom-right of query rectangle
        bottom_right = grid[east - 1][north - 1];

        // Top-right of query rectangle
        if (south - 2 < 0) {
            top_right = 0;
        } else {
            top_right = grid[east - 1][south - 2];
        }

        // Bottom-left of query rectangle
        if (west - 2 < 0) {
            bottom_left = 0;
        } else {
            bottom_left = grid[west - 2][north - 1];
        }

        // Top-left of query rectangle
        if (south - 2 < 0 || west - 2 < 0) {
            top_left = 0;
        } else {
            top_left = grid[west - 2][south - 2];
        }

        // Store results
        result[0] = bottom_right - top_right - bottom_left + top_left; // Population in query rectangle
        result[1] = result[0] / (grid[row_count - 1][column_count - 1]); // Percentage of population in query rectangle

        return result;
    }
}

class V5_CreateGrid extends java.lang.Thread {
    // Member variables
    CensusData census_data;
    int row_count;
    int column_count;
    float[] corners;
    int west;
    int south;
    int east;
    int north;

    int[][] grid;
    int[][] locks;

    int high;
    int low;

    // Constructor
    public V5_CreateGrid(CensusData census_d, int row_c, int col_c, float[] corns, int w, int s, int e, int n,
            int[][] grd, int[][] lcks, int lo, int hi) {
        census_data = census_d;
        row_count = row_c;
        column_count = col_c;
        corners = corns;
        west = w;
        south = s;
        east = e;
        north = n;

        grid = grd;
        locks = lcks;

        high = hi;
        low = lo;
    }

    // Override the run() function
    public void run() {
        for (int i = low; i < high; i++) {
            int blockX = (int) Math.floor((census_data.data[i].longitude - west) / ((east - west) / row_count));
            int blockY = (int) Math.floor((census_data.data[i].latitude - south) / ((north - south) / column_count));
            if (blockX == row_count && blockY == column_count && locks[blockX - 1][blockY - 1] == 0) {
                locks[blockX - 1][blockY - 1] = 1;
                grid[blockX - 1][blockY - 1] += census_data.data[i].population;
                locks[blockX - 1][blockY - 1] = 0;
            } else if (blockX == row_count && locks[blockX - 1][blockY] == 0) {
                locks[blockX - 1][blockY] = 1;
                grid[blockX - 1][blockY] += census_data.data[i].population;
                locks[blockX - 1][blockY] = 0;
            } else if (blockY == column_count && locks[blockX][blockY - 1] == 0) {
                locks[blockX][blockY - 1] = 1;
                grid[blockX][blockY - 1] += census_data.data[i].population;
                locks[blockX][blockY - 1] = 0;
            } else if (locks[blockX][blockY] == 0) {
                locks[blockX][blockY] = 1;
                grid[blockX][blockY] += census_data.data[i].population;
                locks[blockX][blockY] = 1;
            }
        }
    }
}

// //////////////////////////////////////////////////

public class PopulationQuery {
    // next four constants are relevant to parsing
    public static final int TOKENS_PER_LINE = 7;
    public static final int POPULATION_INDEX = 4; // zero-based indices
    public static final int LATITUDE_INDEX = 5;
    public static final int LONGITUDE_INDEX = 6;

    // parse the input file into a large array held in a CensusData object
    public static CensusData parse(String filename) {
        CensusData result = new CensusData();

        try {
            BufferedReader fileIn = new BufferedReader(new FileReader(filename));

            // Skip the first line of the file
            // After that each line has 7 comma-separated numbers (see constants above)
            // We want to skip the first 4, the 5th is the population (an int)
            // and the 6th and 7th are latitude and longitude (floats)
            // If the population is 0, then the line has latitude and longitude of +.,-.
            // which cannot be parsed as floats, so that's a special case
            // (we could fix this, but noisy data is a fact of life, more fun
            // to process the real data as provided by the government)

            String oneLine = fileIn.readLine(); // skip the first line

            // read each subsequent line and add relevant data to a big array
            while ((oneLine = fileIn.readLine()) != null) {
                String[] tokens = oneLine.split(",");
                if (tokens.length != TOKENS_PER_LINE)
                    throw new NumberFormatException();
                int population = Integer.parseInt(tokens[POPULATION_INDEX]);
                if (population != 0)
                    result.add(population, Float.parseFloat(tokens[LATITUDE_INDEX]),
                            Float.parseFloat(tokens[LONGITUDE_INDEX]));
            }

            fileIn.close();
        } catch (IOException ioe) {
            System.err.println("Error opening/reading/writing input or output file.");
            System.exit(1);
        } catch (NumberFormatException nfe) {
            System.err.println(nfe.toString());
            System.err.println("Error in file format");
            System.exit(1);
        }
        return result;
    }

    // ==================================================
    // VERSION 1
    // ==================================================

    // Gets the total population from census data (currently not used)
    private static int v1_GetTotalPopulation(CensusData census_data) {
        int result = 0;
        for (int i = 0; i < census_data.data_size; i++) {
            result += census_data.data[i].population;
        }
        return result;
    }

    // Gets the x, y coordinates of a block
    private static float[] v1_GetBlockCoordinates(CensusData census_data, float[] corners, int row_count,
            int column_count, int index) {
        float result[] = new float[2];
        result[0] = (census_data.data[index].longitude - corners[0]) / ((corners[2] - corners[0]) / row_count) + 1;
        result[1] = (census_data.data[index].latitude - corners[1]) / ((corners[3] - corners[1]) / column_count) + 1;
        return result;
    }

    // Gets the extremes of the lagitude and longitude from the census data
    private static float[] v1_GetCorners(CensusData census_data) {
        // Use the first value of the census data as a starting point
        float north = census_data.data[0].latitude;
        float south = census_data.data[0].latitude;
        float east = census_data.data[0].longitude;
        float west = census_data.data[0].longitude;

        // For each side, assign new value if it is greater or less than the previous
        // value
        for (int i = 0; i < census_data.data_size; i++) {
            float latitude = census_data.data[i].latitude;
            float longitude = census_data.data[i].longitude;
            if (latitude > north)
                north = latitude;
            if (latitude < south)
                south = latitude;
            if (longitude > east)
                east = longitude;
            if (longitude < west)
                west = longitude;
        }

        // After going through all data points, return the results
        float[] result = { west, south, east, north };
        return result;

    }

    // Gets the query results in form of the population within the query rectangle
    // and the percentage of the population within that rectangle compared to total
    // population
    private static float[] v1_GetQueryResults(CensusData census_data, int row_count, int column_count, float[] corners,
            int west, int south, int east, int north) {

        float[] result = new float[2]; // Result to be returned, 1st value is population within query rectanle, 2nd
                                       // value is population percentage
        float total_population = 0; // To count total population

        // Loop over each census data block
        for (int i = 0; i < census_data.data_size; i++) {
            total_population += census_data.data[i].population; // Count the total population

            // Get the x and y values of each census block and add its population if its in
            // the query rectangle
            float[] block_position = v1_GetBlockCoordinates(census_data, corners, row_count, column_count, i);
            int blockX = (int) block_position[0];
            int blockY = (int) block_position[1];

            // System.out.println(blockX + " " + blockY);

            if (blockX >= west && blockX < east + 1 && blockY >= south && blockY < north + 1) {
                result[0] += census_data.data[i].population;
            }
        }

        result[1] = result[0] / total_population; // Calculate and store population percentage
        return result;
    }

    // ==================================================

    // ==================================================
    // VERSION 3
    // ==================================================

    // Proprocesses grid to make queries faster
    private static int[][] v3_PreprocessGrid(CensusData census_data, int row_count, int column_count, float[] corners) {
        int[][] result = v3_CreateGrid(census_data, row_count, column_count, corners);
        result = v3_ProcessGrid(result);
        return result;
    }

    // Fills the grid with the population per block
    public static int[][] v3_CreateGrid(CensusData census_data, int row_count, int column_count, float[] corners) {
        int[][] result = new int[row_count][column_count];

        for (int i = 0; i < census_data.data_size; i++) {
            // Get block coords from each CensusData point
            int blockX = (int) ((census_data.data[i].longitude - corners[0]) / ((corners[2] - corners[0]) / row_count));
            int blockY = (int) ((census_data.data[i].latitude - corners[1])
                    / ((corners[3] - corners[1]) / column_count));

            // Convert block coords to array indices
            int blockX_i = blockX - 1;
            int blockY_i = blockY - 1;

            // Use array indices at edges to avoid out of bounds
            if (blockX == row_count && blockY == column_count) { // Top right corner
                result[blockX_i][blockY_i] += census_data.data[i].population;
            }
            // East edge
            else if (blockX == row_count) {
                result[blockX_i][blockY] += census_data.data[i].population;
            }
            // north edge
            else if (blockY == column_count) {
                result[blockX][blockY_i] += census_data.data[i].population;
            }
            // Not on east or north edge
            else {
                result[blockX][blockY] += census_data.data[i].population;
            }
        }

        // Process the grid after filling it
        v3_ProcessGrid(result);

        // Return the final grid
        return result;
    }

    // Process the grid to hold the total for all positions that are neither farther
    // east or farther south
    public static int[][] v3_ProcessGrid(int[][] grid) {
        // Get the row and column counts
        int row_count = grid.length;
        int column_count = grid[0].length;

        // // Process grid in single pass (doesn't work properly)
        // int[][] orig = grid;
        // for (int i = 1; i < row_count - 1; i++) {
        // for (int j = 1; j < column_count - 1; j++) {
        // grid[i][j] = orig[i][j] + grid[i - 1][j] + grid[i][j + 1] - grid[i - 1][j +
        // 1];
        // }
        // }

        // Iterate over each row starting at the 2nd row and add the previous row
        for (int i = 1; i < row_count; i++) {
            grid[i][0] += grid[i - 1][0];
        }

        // Iterate over each column starting at the 2nd column and add the previous
        // column
        for (int j = 1; j < column_count; j++) {
            grid[0][j] += grid[0][j - 1];
        }

        // Add the previous row, column, and top-left block to current block
        for (int r = 1; r < row_count; r++) {
            for (int c = 1; c < column_count; c++) { // Hehe C++ joke....
                grid[r][c] += grid[r - 1][c] + grid[r][c - 1] - grid[r - 1][c - 1];
            }
        }

        // Return thy grid
        return grid;
    }

    // Gets the reults of the query rectangle
    public static float[] v3_GetQueryResults(int[][] grid, int west, int south, int east, int north, int row_count,
            int column_count) {

        float[] result = new float[2];
        int bottom_right;
        int top_right;
        int bottom_left;
        int top_left;

        // Bottom-right of query rectangle
        bottom_right = grid[east - 1][north - 1];

        // Top-right of query rectangle
        if (south - 2 < 0) {
            top_right = 0;
        } else {
            top_right = grid[east - 1][south - 2];
        }

        // Bottom-left of query rectangle
        if (west - 2 < 0) {
            bottom_left = 0;
        } else {
            bottom_left = grid[west - 2][north - 1];
        }

        // Top-left of query rectangle
        if (south - 2 < 0 || west - 2 < 0) {
            top_left = 0;
        } else {
            top_left = grid[west - 2][south - 2];
        }

        // Store results
        result[0] = bottom_right - top_right - bottom_left + top_left; // Population in query rectangle
        result[1] = result[0] / (grid[row_count - 1][column_count - 1]); // Percentage of population in query rectangle

        return result;
    }

    // ==================================================
    // MAIN FUNCTION
    // ==================================================

    // argument 1: file name for input data: pass this to parse
    // argument 2: number of x-dimension buckets
    // argument 3: number of y-dimension buckets
    // argument 4: -v1, -v2, -v3, -v4, or -v5
    public static void main(String[] args) {
        long time_start;
        long time_end;

        // Read and store the census data
        CensusData census_data = new CensusData();
        census_data = parse(args[0]);

        // Store number of total rows and columns of the grid
        int row_count = Integer.parseInt(args[1]);
        int column_count = Integer.parseInt(args[2]);

        ///////////////
        // Version 1
        if (args[3].equals("-v1")) {
            float[] corners = v1_GetCorners(census_data); // Get corners in order of: {west, south, east, north}
            // System.out.println(corners[0] + " " + corners[1] + " " + corners[2] + " " +
            // corners[3]);

            boolean exit = false;
            while (!exit) {
                // Input query
                System.out.println("Please give west, south, east, north coordinates of your query rectangle:");
                String query = System.console().readLine();

                // Split the values using whitespaces, convert them to ints, and store them
                // separately
                String[] query_split = query.split("\\s+");
                int west = Integer.parseInt(query_split[0]);
                int south = Integer.parseInt(query_split[1]);
                int east = Integer.parseInt(query_split[2]);
                int north = Integer.parseInt(query_split[3]);

                // Exit condition
                if (query_split.length != 4) {
                    exit = true;
                    System.exit(1);
                }

                float[] query_results = v1_GetQueryResults(census_data, row_count, column_count, corners, west, south,
                        east, north); // Get results from the query rectangle in order of: {population, percentage}

                System.out.printf("population of rectangle: %d\npercent of total population: %.2f\n",
                        (int) query_results[0], query_results[1] * 100);
            }
        }

        ///////////////
        // Version 2
        else if (args[3].equals("-v2")) {
            float[] corners = V2_GetCorners.getCorners(census_data); // Get corners in order of: {west, south, east,
                                                                     // north}
            // System.out.println(corners[0] + " " + corners[1] + " " + corners[2] + " " +
            // corners[3]);

            boolean exit = false;
            while (!exit) {
                // Input query
                System.out.println("Please give west, south, east, north coordinates of your query rectangle:");
                String query = System.console().readLine();

                // Split the values using whitespaces, convert them to ints, and store them
                // separately
                String[] query_split = query.split("\\s+");
                int west = Integer.parseInt(query_split[0]);
                int south = Integer.parseInt(query_split[1]);
                int east = Integer.parseInt(query_split[2]);
                int north = Integer.parseInt(query_split[3]);

                // Exit condition
                if (query_split.length != 4) {
                    exit = true;
                    System.exit(1);
                }

                float[] query_results = V2_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                        corners, west, south, east, north); // Get results from the query rectangle in order of:
                                                            // {population,
                                                            // percentage}

                System.out.printf("population of rectangle: %d\npercent of total population: %.2f\n",
                        (int) query_results[0], query_results[1] * 100);
            }
        }

        ///////////////
        // Version 3
        else if (args[3].equals("-v3")) {
            float[] corners = v1_GetCorners(census_data); // Get corners in order of: {west, south, east, north}
            // System.out.println(corners[0] + " " + corners[1] + " " + corners[2] + " " +
            // corners[3]);

            // Preprocess data
            // int[][] grid = v3_PreprocessGrid(census_data, row_count, column_count,
            // corners);
            int[][] grid = v3_CreateGrid(census_data, row_count, column_count, corners);

            boolean exit = false;
            while (!exit) {
                // Input query
                System.out.println("Please give west, south, east, north coordinates of your query rectangle:");
                String query = System.console().readLine();

                // Split the values using whitespaces, convert them to ints, and store them
                // separately
                String[] query_split = query.split("\\s+");
                int west = Integer.parseInt(query_split[0]);
                int south = Integer.parseInt(query_split[1]);
                int east = Integer.parseInt(query_split[2]);
                int north = Integer.parseInt(query_split[3]);

                // Exit condition
                if (query_split.length != 4) {
                    exit = true;
                    System.exit(1);
                }

                float[] query_results = v3_GetQueryResults(grid, west, south, east, north, row_count, column_count);

                System.out.printf("population of rectangle: %d\npercent of total population: %.2f\n",
                        (int) query_results[0], query_results[1] * 100);
            }
        }

        ///////////////
        // Version 4
        else if (args[3].equals("-v4")) {
            float[] corners = V2_GetCorners.getCorners(census_data); // Get corners

            boolean exit = false;
            while (!exit) {
                // Input query
                System.out.println("Please give west, south, east, north coordinates of your query rectangle:");
                String query = System.console().readLine();

                // Split the values using whitespaces, convert them to ints, and store them
                String[] query_split = query.split("\\s+");
                int west = Integer.parseInt(query_split[0]);
                int south = Integer.parseInt(query_split[1]);
                int east = Integer.parseInt(query_split[2]);
                int north = Integer.parseInt(query_split[3]);

                // Exit condition
                if (query_split.length != 4) {
                    exit = true;
                    System.exit(1);
                }

                float[] query_results = V4_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                        corners, west, south, east, north);

                System.out.printf("population of rectangle: %d\npercent of total population: %.2f\n",
                        (int) query_results[0], query_results[1] * 100);
            }
        }

        ///////////////
        // Version 5
        else if (args[3].equals("-v5")) {
            float[] corners = V2_GetCorners.getCorners(census_data); // Get corners

            boolean exit = false;
            while (!exit) {
                // Input query
                System.out.println("Please give west, south, east, north coordinates of your query rectangle:");
                String query = System.console().readLine();

                // Split the values using whitespaces, convert them to ints, and store them
                String[] query_split = query.split("\\s+");
                int west = Integer.parseInt(query_split[0]);
                int south = Integer.parseInt(query_split[1]);
                int east = Integer.parseInt(query_split[2]);
                int north = Integer.parseInt(query_split[3]);

                // Exit condition
                if (query_split.length != 4) {
                    exit = true;
                    System.exit(1);
                }

                // Input query
                System.out.println("How many threads would you like to use?:");
                String thread_count_response = System.console().readLine();
                int thread_count = Integer.parseInt(thread_count_response);

                V5_GetQueryResults V5_queryResults = new V5_GetQueryResults(census_data, row_count, column_count,
                        corners, west, south, east, north, thread_count);

                float[] query_results = V5_queryResults.V5_getQueryResults();

                System.out.printf("population of rectangle: %d\npercent of total population: %.2f\n",
                        (int) query_results[0], query_results[1] * 100);
            }
        }

        ///////////////
        // Testing
        else if (args[3].equals("-t")) {

            String query = "1 1 1000 1000";

            ///////////////
            // Version 1
            {
                time_start = System.currentTimeMillis();
                float[] corners = v1_GetCorners(census_data);
                time_end = System.currentTimeMillis();
                System.out.println("v1 corners (ms): " + (time_end - time_start));

                // Trial 1
                {
                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v1_GetQueryResults(census_data, row_count, column_count, corners, west,
                            south, east, north); // Get results from the query rectangle in order of: {population,
                                                 // percentage}
                    time_end = System.currentTimeMillis();
                    System.out.println("v1 query results (ms): " + (time_end - time_start));
                }

                // Trial 2
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v1_GetQueryResults(census_data, row_count, column_count, corners, west,
                            south, east, north); // Get results from the query rectangle in order of: {population,
                                                 // percentage}
                    time_end = System.currentTimeMillis();
                    System.out.println("v1 query results (ms): " + (time_end - time_start));
                }

                // Trial 3
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v1_GetQueryResults(census_data, row_count, column_count, corners, west,
                            south, east, north); // Get results from the query rectangle in order of: {population,
                                                 // percentage}
                    time_end = System.currentTimeMillis();
                    System.out.println("v1 query results (ms): " + (time_end - time_start));
                }

                // Tiral 4
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v1_GetQueryResults(census_data, row_count, column_count, corners, west,
                            south, east, north); // Get results from the query rectangle in order of: {population,
                                                 // percentage}
                    time_end = System.currentTimeMillis();
                    System.out.println("v1 query results (ms): " + (time_end - time_start));
                }

                // Trial 5
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v1_GetQueryResults(census_data, row_count, column_count, corners, west,
                            south, east, north); // Get results from the query rectangle in order of: {population,
                                                 // percentage}
                    time_end = System.currentTimeMillis();
                    System.out.println("v1 query results (ms): " + (time_end - time_start));
                }
            }

            ///////////////
            // Version 2
            {
                time_start = System.currentTimeMillis();
                float[] corners = v1_GetCorners(census_data);
                time_end = System.currentTimeMillis();
                System.out.println("v2 corners (ms): " + (time_end - time_start));

                // Trial 1
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V2_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v2 query results (ms): " + (time_end - time_start));
                }

                // Trial 2
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V2_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v2 query results (ms): " + (time_end - time_start));
                }

                // Trial 3
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V2_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v2 query results (ms): " + (time_end - time_start));
                }

                // Tiral 4
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V2_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v2 query results (ms): " + (time_end - time_start));
                }

                // Trial 5
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V2_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v2 query results (ms): " + (time_end - time_start));
                }
            }
            ///////////////
            // Version 3
            {
                time_start = System.currentTimeMillis();
                float[] corners = v1_GetCorners(census_data);
                time_end = System.currentTimeMillis();
                System.out.println("v3 corners (ms): " + (time_end - time_start));

                time_start = System.currentTimeMillis();
                int[][] grid = v3_CreateGrid(census_data, row_count, column_count, corners);
                time_end = System.currentTimeMillis();
                System.out.println("v3 create grid (ms): " + (time_end - time_start));

                // Trial 1
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v3_GetQueryResults(grid, west, south, east, north, row_count, column_count);
                    time_end = System.currentTimeMillis();
                    System.out.println("v3 query results (ms): " + (time_end - time_start));
                }

                // Trial 2
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v3_GetQueryResults(grid, west, south, east, north, row_count, column_count);
                    time_end = System.currentTimeMillis();
                    System.out.println("v3 query results (ms): " + (time_end - time_start));
                }

                // Trial 3
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v3_GetQueryResults(grid, west, south, east, north, row_count, column_count);
                    time_end = System.currentTimeMillis();
                    System.out.println("v3 query results (ms): " + (time_end - time_start));
                }

                // Tiral 4
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v3_GetQueryResults(grid, west, south, east, north, row_count, column_count);
                    time_end = System.currentTimeMillis();
                    System.out.println("v3 query results (ms): " + (time_end - time_start));
                }

                // Trial 5
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = v3_GetQueryResults(grid, west, south, east, north, row_count, column_count);
                    time_end = System.currentTimeMillis();
                    System.out.println("v3 query results (ms): " + (time_end - time_start));
                }
            }

            ///////////////
            // Version 4
            {
                time_start = System.currentTimeMillis();
                float[] corners = V2_GetCorners.getCorners(census_data); // Get corners
                time_end = System.currentTimeMillis();
                System.out.println("v4 corners (ms): " + (time_end - time_start));

                // Trial 1
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V4_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v4 query results (ms): " + (time_end - time_start));

                }

                // Trial 2
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V4_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v4 query results (ms): " + (time_end - time_start));
                }

                // Trial 3
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V4_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v4 query results (ms): " + (time_end - time_start));
                }

                // Tiral 4
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V4_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v4 query results (ms): " + (time_end - time_start));
                }

                // Trial 5
                {

                    String[] query_split = query.split("\\s+");
                    int west = Integer.parseInt(query_split[0]);
                    int south = Integer.parseInt(query_split[1]);
                    int east = Integer.parseInt(query_split[2]);
                    int north = Integer.parseInt(query_split[3]);

                    time_start = System.currentTimeMillis();
                    float[] query_results = V4_GetQueryResults.getQueryResults(census_data, row_count, column_count,
                            corners, west, south, east, north);
                    time_end = System.currentTimeMillis();
                    System.out.println("v4 query results (ms): " + (time_end - time_start));
                }
            }

        }

        else {
            System.out.println("Argument 4 invalid");
            System.exit(1);
        }

    }
}
