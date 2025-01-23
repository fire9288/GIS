#pragma warning(disable:4996)
/**
 * This program implements a landscape imagery analysis system with data compression
 * capabilities. It processes geographic data using three different compression methods
 * and provides functionality for querying landscape visibility and locations.
 */

#include <iostream>
#include <stdio.h>
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_conv.h"
#include<math.h>
#include <stdlib.h>
#include <time.h>
#include<Windows.h>
#include<tchar.h>
#include<Psapi.h>
using namespace std;

#define NPoint 70 // Number of observation points/layers in the landscape data

typedef unsigned char U_CHAR;

typedef struct Coordinate {
    double x;  // X coordinate (longitude/easting)
    double y;  // Y coordinate (latitude/northing)
} Coordinate;

// Structure for storing grid coordinates in a linked list
typedef struct coord {
    int r, c;           // Row and column indices
    struct coord* next; // Pointer to next coordinate
} coord;

// Structure for storing compressed landscape data
typedef struct Node {
    unsigned int* arr;   // Array storing compressed data
    struct Node* down;   // Pointer to next node
    coord* right;        // Pointer to coordinate list
} Node;

// Structure for managing compressed area data
typedef struct Area {
    int len;           // Length of compressed data array
    Node* point;       // Pointer to first node
    int size;         // Number of unique patterns
} Area;

/**
 * Compares two unsigned int arrays for equality
 * @param a First array
 * @param b Second array
 * @param len Length of arrays
 * @return true if arrays are equal, false otherwise
 */
bool compare(unsigned int* a, unsigned int* b, int len);

/**
 * Searches for a matching data pattern in the compressed area structure
 * @param ar Pointer to Area structure
 * @param data Data pattern to search for
 * @return Pointer to matching Node if found, NULL otherwise
 */
Node* lookfor(Area* ar, unsigned int* data);

/**
 * First compression stage: Compresses 3D binary data into bit-packed format
 * @param x Input 3D array
 * @param len Number of layers
 * @param M Number of rows
 * @param N Number of columns
 * @return Compressed 3D array
 */
unsigned int*** compress1(unsigned int*** x, int len, int M, int N);

/**
 * Writes first compression stage results to file
 */
void write1(unsigned int*** goal, int len, int M, int N);

/**
 * Identifies unique patterns in compressed data
 * @return Area structure containing unique patterns
 */
Area* unique(unsigned int*** data, int len, int M, int N);

/**
 * Second compression stage: Creates lookup table for unique patterns
 */
unsigned int** compress2(Area* comArea, unsigned int** data);

/**
 * Writes second compression stage results to file
 */
void write2(Area* comArea, unsigned int** data_status, unsigned int** data, int M, int N);

/**
 * Third compression stage: Implements run-length encoding
 */
void write3(Area* comArea, unsigned int** data_status, unsigned int** a, int M, int N);


/**
 * Reads the original uncompressed data from benchmark.txt file
 * @param goal Reference to 3D array where data will be stored
 * @param len Number of layers in the data
 * @param M Number of rows
 * @param N Number of columns
 */
void read0(unsigned int***& goal, int len, int M, int N);

/**
 * Reads the first compression stage data from firstCompress.txt file
 * This contains the bit-packed representation of the original data
 * @param goal Reference to 3D array where compressed data will be stored
 * @param len Number of layers after bit-packing (len = original_len/32 + remainder)
 * @param M Number of rows
 * @param N Number of columns
 */
void read1(unsigned int***& goal, int len, int M, int N);

/**
 * Reads the second compression stage data from secondCompress.txt file
 * This includes both the lookup table of unique patterns and the grid references
 * @param comArea Reference to Area structure containing pattern information
 * @param data_status Reference to 2D array storing the lookup table
 * @param data Reference to 2D array storing grid references
 * @param M Number of rows
 * @param N Number of columns
 */
void read2(Area*& comArea, unsigned int**& data_status, unsigned int**& data, int M, int N);

/**
 * Reads the third compression stage data from thirdCompress.txt file
 * This includes the run-length encoded representation of the grid data
 * @param comArea Reference to Area structure containing pattern information
 * @param data_status Reference to 2D array storing the lookup table
 * @param a Reference to 2D array where decoded run-length data will be stored
 * @param M Number of rows
 * @param N Number of columns
 */
void read3(Area*& comArea, unsigned int**& data_status, unsigned int**& a, int M, int N);

/**
 * Compares two arrays of unsigned integers for equality
 * @param a First array to compare
 * @param b Second array to compare
 * @param len Length of both arrays
 * @return 1 if arrays are equal, 0 if they differ or if either is NULL
 *
 * Example usage:
 * unsigned int arr1[] = {1, 2, 3};
 * unsigned int arr2[] = {1, 2, 3};
 * if (equals(arr1, arr2, 3)) {
 *     // Arrays are equal
 * }
 */
int equals(unsigned int* a, unsigned int* b, int len)
{
    int i = 0;
    if (NULL == a || NULL == b)
        return 0;
    while (i < len)
    {
        if (a[i] != b[i])
            return 0;
        i++;
    }
    return 1;
}

/**
 * Decodes a bit-packed integer array back to binary array
 * @param x Input packed array
 * @param n Length of output array
 * @return Decoded binary array
 */
int* decode(unsigned int x[], int n)
{
    int i, k;
    int* ch;
    ch = (int*)malloc(n * sizeof(int));
    for (i = 0; i < n; i++)
        ch[i] = 0;
    k = 0;
    int lencomp = n / 32 + (n % 32 == 0 ? 0 : 1);
    for (i = 0; i < lencomp; i++)
    {
        unsigned int temp = x[i];
        k = i * 32;
        while (temp != 0)
        {
            ch[k] = (temp % 2 == 1) ? 1 : 0;
            temp = temp / 2;
            k++;
        }
    }
    return ch;
}

/**
 * Encodes a binary array into bit-packed format
 * @param temp Input binary array
 * @param len Length of input array
 * @return Packed integer array
 */
unsigned int* encode(int* temp, int len)//Compress encoding into integer
{
    //cout << "hello\n" << endl;
    int lenComp = len / 32 + (len % 32 == 0 ? 0 : 1);
    unsigned int* tempIntcode = (unsigned int*)malloc(lenComp * sizeof(unsigned int));
    unsigned int t = 0;
    int k = 0;
    for (int i = 0; i < lenComp && k < len; i++)
    {
        t = 0;
        for (int j = 0; j < 32; j++)
        {
            if (temp[k] == 1)
            {
                int js = 1;
                for (int r = 1; r <= j; r++)
                {
                    js *= 2;
                }
                t += js;
            }
            k++;
        }
        tempIntcode[i] = t;
    }
    return tempIntcode;
}

/**
 * Queries landscape imagery visibility for original uncompressed data
 */
int* inquire_landscapeImagery0(unsigned int*** x, int len, int M, int N, Coordinate left, Coordinate right, int coorMinX, int coorMinY, double resolution)
{
    int rowStart, rowEnd, colStart, colEnd;
    double minx, miny, maxx, maxy;
    double leftb = coorMinX, lowb = coorMinY;
    int* temp = (int*)malloc(len * sizeof(int));
    for (int i = 0; i < len; i++)
        temp[i] = 0;
    if (left.x < right.x)
    {
        minx = left.x;
        maxx = right.x;
    }
    else
    {
        maxx = left.x;
        minx = right.x;
    }
    if (left.y < right.y)
    {
        miny = left.y;
        maxy = right.y;
    }
    else
    {
        maxy = left.y;
        miny = right.y;
    }
    rowStart = (M - 1) - (int)((maxy - lowb) / resolution);
    rowEnd = (M - 1) - (int)((miny - lowb) / resolution);
    colStart = (int)((minx - leftb) / resolution);
    colEnd = (int)((maxx - leftb) / resolution);
    for (int i = rowStart; i < rowEnd; i++)
    {
        for (int j = colStart; j < colEnd; j++)
        {

            for (int k = 0; k < len; k++)
            {
                temp[k] = temp[k] | x[k][i][j];
            }

        }
    }
    return temp;
}

/**
 * Queries landscape imagery visibility for first compression stage
 */
int* inquire_landscapeImagery1(unsigned int*** x, int lencomp, int len, int M, int N, Coordinate left, Coordinate right, int coorMinX, int coorMinY, double resolution)
{
    int rowStart, rowEnd, colStart, colEnd;
    double minx, miny, maxx, maxy;
    double leftb = coorMinX, lowb = coorMinY;
    unsigned int* temp = (unsigned int*)malloc(lencomp * sizeof(unsigned int));
    for (int i = 0; i < lencomp; i++)
        temp[i] = 0;
    if (left.x < right.x)
    {
        minx = left.x;
        maxx = right.x;
    }
    else
    {
        maxx = left.x;
        minx = right.x;
    }
    if (left.y < right.y)
    {
        miny = left.y;
        maxy = right.y;
    }
    else
    {
        maxy = left.y;
        miny = right.y;
    }
    rowStart = (M - 1) - (int)((maxy - lowb) / resolution);
    rowEnd = (M - 1) - (int)((miny - lowb) / resolution);
    colStart = (int)((minx - leftb) / resolution);
    colEnd = (int)((maxx - leftb) / resolution);
    for (int i = rowStart; i < rowEnd; i++)
    {
        for (int j = colStart; j < colEnd; j++)
        {

            for (int k = 0; k < lencomp; k++)
            {
                temp[k] = temp[k] | x[k][i][j];
            }

        }
    }
    int* p = decode(temp, len);
    return p;
}

/**
 * Queries landscape imagery visibility for second compression stage
 */
int* inquire_landscapeImagery2(unsigned int** x, unsigned int indexStatus, unsigned int** dataStatus, int lencomp, int len, int M, int N, Coordinate left, Coordinate right, int coorMinX, int coorMinY, double resolution)
{
    int rowStart, rowEnd, colStart, colEnd;
    double minx, miny, maxx, maxy;
    double leftb = coorMinX, lowb = coorMinY;
    unsigned int* temp = (unsigned int*)malloc(lencomp * sizeof(unsigned int));
    for (int i = 0; i < lencomp; i++)
        temp[i] = 0;
    if (left.x < right.x)
    {
        minx = left.x;
        maxx = right.x;
    }
    else
    {
        maxx = left.x;
        minx = right.x;
    }
    if (left.y < right.y)
    {
        miny = left.y;
        maxy = right.y;
    }
    else
    {
        maxy = left.y;
        miny = right.y;
    }
    rowStart = (M - 1) - (int)((maxy - lowb) / resolution);
    rowEnd = (M - 1) - (int)((miny - lowb) / resolution);
    colStart = (int)((minx - leftb) / resolution);
    colEnd = (int)((maxx - leftb) / resolution);
    for (int i = rowStart; i < rowEnd; i++)
    {
        for (int j = colStart; j < colEnd; j++)
        {

            for (int k = 0; k < lencomp; k++)
            {
                temp[k] = temp[k] | dataStatus[x[i][j]][k];
            }

        }
    }
    int* p = decode(temp, len);
    return p;
}

/**
 * Queries locations matching specified landscape imagery pattern in original uncompressed data
 *
 * @param x 3D array of original landscape data [layer][row][col]
 * @param len Number of layers in the data
 * @param M Number of rows
 * @param N Number of columns
 * @param landscapeImagery Binary array specifying which layers to consider (1 = consider, 0 = ignore)
 * @return 2D array where 1 indicates matching locations and 0 indicates non-matching locations
 *
 * Algorithm:
 * 1. Initialize result array with all 1's
 * 2. For each layer marked as 1 in landscapeImagery:
 *    - Update result using logical AND with that layer's data
 * 3. Final result shows locations matching all selected layers
 */
int** inquire_location0(unsigned int*** x, int len, int M, int N, int* landscapeImagery)
{
    int** p = (int**)malloc(M * sizeof(int*));
    for (int i = 0; i < M; i++)
    {
        p[i] = (int*)malloc(N * sizeof(int));
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
            p[i][j] = 1;
    }
    for (int k = 0; k < len; k++)
    {
        if (landscapeImagery[k] == 1)
        {
            for (int i = 0; i < M; i++)
                for (int j = 0; j < N; j++)
                    p[i][j] = p[i][j] && x[k][i][j];
        }
    }
    return p;
}

/**
 * Queries locations matching specified landscape imagery pattern in first compression stage data
 *
 * @param x 3D array of compressed landscape data
 * @param lencomp Length of compressed data (original_length/32 + remainder)
 * @param M Number of rows
 * @param N Number of columns
 * @param landscapeImagery Compressed pattern to search for
 * @return 2D array where 1 indicates matching locations and 0 indicates non-matching locations
 *
 * Algorithm:
 * 1. Initialize result array with all 0's
 * 2. For each grid position:
 *    - Compare compressed data with target pattern
 *    - Mark matching positions with 1
 */
int** inquire_location1(unsigned int*** x, int lencomp, int M, int N, unsigned int* landscapeImagery)
{
    int** p = (int**)malloc(M * sizeof(int*));
    unsigned int* temp = (unsigned int*)malloc(sizeof(unsigned int) * lencomp);
    for (int i = 0; i < M; i++)
    {
        p[i] = (int*)malloc(N * sizeof(int));
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
            p[i][j] = 0;
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < lencomp; k++)
                temp[k] = x[k][i][j];
            if (equals(temp, landscapeImagery, lencomp))
                p[i][j] = 1;
        }
    }
    return p;
}

/**
 * Queries locations matching specified index in second compression stage data
 *
 * @param x 2D array of indices referencing the compressed patterns
 * @param M Number of rows
 * @param N Number of columns
 * @param indexLandscapeImagery Index of the pattern to search for
 * @return 2D array where 1 indicates matching locations and 0 indicates non-matching locations
 *
 * Algorithm:
 * 1. Initialize result array with all 0's
 * 2. Mark all positions matching the target index with 1
 *
 * Note: There appears to be a bug in the comparison - using = instead of ==
 * Current: if (x[i][j] = indexLandscapeImagery)
 * Should be: if (x[i][j] == indexLandscapeImagery)
 */
int** inquire_location2(unsigned int** x, int M, int N, unsigned int indexLandscapeImagery)
{
    int** p = (int**)malloc(M * sizeof(int*));
    for (int i = 0; i < M; i++)
    {
        p[i] = (int*)malloc(N * sizeof(int));
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
            p[i][j] = 0;
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (x[i][j] = indexLandscapeImagery)
                p[i][j] = 1;
        }
    }
    return p;
}

/**
 * Main function to demonstrate landscape data compression and analysis
 * Processes multiple GeoTIFF files, applies three compression stages,
 * and performs performance testing of visibility and location queries
 */
int main() {
    clock_t start, finish;
    double duration;
    FILE* file;

    GDALAllRegister();// Initialize GDAL for GeoTIFF processing
    GDALDataset* poDataset;
    double M_x = 513249, M_y = 4.50014e+06;  // Study area minimum coordinates
    int MM = 840, NN = 1200;       // Height and width of raster

    unsigned int*** arr = (unsigned int***)malloc(NPoint * sizeof(unsigned int**));// Allocate 3D array for multiple viewshed
    unsigned int*** goal;//Compressed 3D array
    int lencompress;//Length of Compressed 3D array
    Area* comArea;//Deduplicated visibility data structure
    unsigned int** data_status;//Visible observation point sequence
    unsigned int** data2;//Correspondence relationship between visible observation point status and location grid

    // Build the necessary data structures
    for (int i = 0; i < NPoint; i++) {
        arr[i] = (unsigned int**)malloc(MM * sizeof(unsigned int*));
        for (int j = 0; j < MM; j++) {
            arr[i][j] = (unsigned int*)malloc(NN * sizeof(unsigned int));
        }
    }
    lencompress = NPoint / 32 + (NPoint % 32 == 0 ? 0 : 1);
    goal = (unsigned int***)malloc(lencompress * sizeof(unsigned int**));
    for (int i = 0; i < lencompress; i++)
    {
        goal[i] = (unsigned int**)malloc(MM * sizeof(unsigned int*));
        for (int j = 0; j < MM; j++)
            goal[i][j] = (unsigned int*)malloc(NN * sizeof(unsigned int));
    }

    // Process each GeoTIFF file
    for (int i = 0; i < NPoint; i++) {
        // Construct input filename
        char strin[5];
        sprintf(strin, "%d", i);
        char file_in[60] = "F:/secondStepOut/viewshed100/viewshed_";
        strcat(file_in, strin);
        strcat(file_in, "_masked.tif");

        // Open and validate GeoTIFF file
        poDataset = static_cast<GDALDataset*>(GDALOpen(file_in, GA_ReadOnly));
        if (poDataset == NULL) {
            std::cout << "Cannot open TIFF file." << std::endl;
            return 1;
        }

        // Get raster dimensions
        int nWidth = poDataset->GetRasterXSize();
        int nHeight = poDataset->GetRasterYSize();
        MM = nHeight;
        NN = nWidth;
        // Allocate memory for current layer
        arr[i] = (unsigned int**)malloc(nHeight * sizeof(unsigned int*));
        for (int j = 0; j < nHeight; j++) {
            arr[i][j] = (unsigned int*)malloc(nWidth * sizeof(unsigned int));
        }

        // Get study area coordinates
        double adfGeoTransform[6];
        poDataset->GetGeoTransform(adfGeoTransform);
        double minX = adfGeoTransform[0];
        double maxY = adfGeoTransform[3];
        double maxX = adfGeoTransform[0] + adfGeoTransform[1] * nWidth + adfGeoTransform[2] * nHeight;
        double minY = adfGeoTransform[3] + adfGeoTransform[4] * nWidth + adfGeoTransform[5] * nHeight;
        M_x = minX;
        M_y = minY;
        // Read first band of TIFF file
        GDALRasterBand* poBand = poDataset->GetRasterBand(1);

        // Process each pixel
        for (int j = 0; j < nHeight; j++) {
            for (int k = 0; k < nWidth; k++) {
                float pixel;
                GDALRasterIO(poBand, GF_Read, k, j, 1, 1, &pixel, 1, 1, GDT_Float32, 0, 0);
                int temp = (int)pixel;
                arr[i][j][k] = temp;
                // Binarize: 1 remains 1, everything else becomes 0
                if (1 != arr[i][j][k])
                    arr[i][j][k] = 0;
            }
        }
        GDALClose(poDataset);
    }

    // Write original 3D data to benchmark file
    file = fopen("f:/outVisibilityData/benchmark.txt", "w");
    for (int i = 0; i < NPoint; i++) {
        for (int j = 0; j < MM; j++) {
            for (int k = 0; k < NN; k++) {
                fprintf(file, "%u ", arr[i][j][k]);
            }
            fprintf(file, "\n\n");
        }
    }
    fclose(file);

    // First compression stage: Bit packing

    start = clock();
    goal = compress1(arr, NPoint, MM, NN);
    lencompress = NPoint / 32 + (NPoint % 32 == 0 ? 0 : 1);
    write1(goal, lencompress, MM, NN);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << endl << "First compression time:" << duration << endl;

    // Second compression stage: Pattern deduplication
    start = clock();
    comArea = unique(goal, lencompress, MM, NN);
    data2 = (unsigned int**)malloc(MM * sizeof(unsigned int*));
    for (int i = 0; i < MM; i++)
        data2[i] = (unsigned int*)malloc(NN * sizeof(unsigned int));
    data_status = compress2(comArea, data2);
    write2(comArea, data_status, data2, MM, NN);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << endl << "Second compression time:" << duration << endl;

    // Third compression stage: Run-length encoding
    start = clock();
    write3(comArea, data_status, data2, MM, NN);
    finish = clock();
    duration = (double)(finish - start) / CLOCKS_PER_SEC;
    cout << endl << "Third compression time:" << duration << endl;

    //// Define test coordinates for visibility analysis
    //Coordinate first[] = { {518800,4506041},{518300,4504041},{516800,4503188},
    //                      {515800,4502400},{517240,4501000} };
    //Coordinate second[] = { {521300,4506666},{521300,4506666},{521300,4506666},
    //                       {521300,4506666},{521210,4508888} };
    //

    //cout << "==Now measuring visibility query performance==" << endl;
    //for (int kk = 0; kk < 5; kk++) {
    //    // Test original data structure
    //    start = clock();
    //    read0(arr, NPoint, MM, NN);        
    //    for (int k = 0; k < 1000; k++)
    //        inquire_landscapeImagery0(arr, NPoint, MM, NN, first[kk], second[kk], M_x, M_y, 12.5);
    //    finish = clock();
    //    duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    //    cout << endl << "benchmark:" << duration << "   ";

    //    // Test first compression
    //    
    //    start = clock();
    //    read1(goal, lencompress, MM, NN);        
    //    for (int k = 0; k < 1000; k++)
    //        inquire_landscapeImagery1(goal, lencompress, NPoint, MM, NN,
    //            first[kk], second[kk], M_x, M_y, 12.5);
    //    finish = clock();
    //    duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    //    cout << endl << "firstCompress:" << duration << "   ";

    //    comArea = unique(goal, lencompress, MM, NN);    
    //    data2 = (unsigned int**)malloc(MM * sizeof(unsigned int*));
    //    for (int i = 0; i < MM; i++)
    //        data2[i] = (unsigned int*)malloc(NN * sizeof(unsigned int));
    //    data_status = compress2(comArea, data2);

    //    start = clock();
    //    read2(comArea, data_status, data2, MM, NN);
    //    for (int k = 0; k < 1000; k++)
    //        inquire_landscapeImagery2(data2, 2, data_status, lencompress, NPoint, MM, NN, first[kk], second[kk], M_x, M_y, 12.5);
    //    /*int* t2 = inquire_landscapeImagery2(data2, 2, data_status, lencompress, NPoint, MM, NN, fir, sec, M_x, M_y, 12.5);
    //    for (int i = 0; i < lencompress; i++)
    //        cout << " " << t2[i];*/
    //    finish = clock();
    //    duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    //    cout << endl << "secondCompress:" << duration;

    //    start = clock();
    //    read3(comArea, data_status, data2, MM, NN);
    //    for (int k = 0; k < 1000; k++)
    //        inquire_landscapeImagery2(data2, 2, data_status, lencompress, NPoint, MM, NN, first[kk], second[kk], M_x, M_y, 12.5);
    //    finish = clock();
    //    duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    //    cout << endl << "thirdCompress:" << duration << endl;
    //}
    //
    //cout << "==Now measuring query location performance==" << endl;
    //int** landscapeIm = (int**)malloc(5 * sizeof(int*));
    //int number = 70;
    //
    //for (int j = 0; j < 5; j++)
    //{
    //    landscapeIm[j] = (int*)malloc(number * sizeof(int));
    //    for (int k = 0; k < number; k++)
    //        landscapeIm[j][k] = 0;
    //    for (int k = 0; k < (j + 1) * 10; k++)
    //        landscapeIm[j][k] = 1;

    //    start = clock();
    //    read0(arr, NPoint, MM, NN);
    //    for (int k = 0; k < 1000; k++)
    //    {
    //        int** temploc = inquire_location0(arr, number, MM, NN, landscapeIm[j]);
    //        for (int q = 0; q < MM; q++)
    //            delete (temploc[q]);
    //        delete temploc;
    //    }
    //    finish = clock();
    //    duration = (double)(finish - start);
    //    cout << endl << "benchmark:" << duration << "   ";
    //    unsigned* landcom = NULL;

    //    if (NULL != landcom)
    //        free(landcom);

    //    landcom = encode(landscapeIm[j], number);
    //    start = clock();
    //    read1(goal, lencompress, MM, NN);
    //    for (int k = 0; k < 1000; k++)
    //    {
    //        int** temploc = inquire_location1(goal, lencompress, MM, NN, landcom);
    //        for (int q = 0; q < MM; q++)
    //            delete (temploc[q]);
    //        delete temploc;
    //    }

    //    finish = clock();
    //    duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    //    cout << endl << "firstCompress:" << duration << "   ";

    //    srand(time(NULL)); 
    //    int randomNumber = (rand() % (10 - 1 + 1)) + 1;         
    //    start = clock();
    //    read2(comArea, data_status, data2, MM, NN);
    //    for (int k = 0; k < 1000; k++)
    //    {
    //        int** temploc = inquire_location2(data2, MM, NN, randomNumber);
    //        for (int q = 0; q < MM; q++)
    //            delete (temploc[q]);
    //        delete temploc;
    //    }
    //    finish = clock();
    //    duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    //    cout << endl << "secondCompress:" << duration ;

    //    start = clock();
    //    read3(comArea, data_status, data2, MM, NN);
    //    for (int k = 0; k < 1000; k++)
    //    {
    //        int** temploc = inquire_location2(data2, MM, NN, randomNumber);
    //        for (int q = 0; q < MM; q++)
    //            delete (temploc[q]);
    //        delete temploc;
    //    }
    //    finish = clock();
    //    duration = (double)(finish - start) / CLOCKS_PER_SEC * 1000;
    //    cout << endl << "thirdCompress:" << duration << endl;
    //}


    //Coordinate Afir = { 519628,4505179 };
    //Coordinate Asec = { 520228,4505599 };
    //
    //Coordinate Bfir = { 520047,4503335 };
    //Coordinate Bsec = { 520647,4503755 };
    //read1(goal, lencompress, MM, NN);
    //int* pA = inquire_landscapeImagery1(goal, lencompress, NPoint, MM, NN, Afir, Asec, M_x, M_y, 12.5);
    //cout << "area A:" << endl;
    //for (int i = 0; i < NPoint; i++)
    //    cout << " " << pA[i];

    //pA = inquire_landscapeImagery1(goal, lencompress, NPoint, MM, NN, Bfir, Bsec, M_x, M_y, 12.5);
    //cout << endl << "area B:" << endl;
    //for (int i = 0; i < NPoint; i++)
    //    cout << " " << pA[i];

    //int landscape_wall[NPoint] = { 0 };
    //int landscape_Tower[NPoint] = { 0 };
    //int Enemy[7] = { 6,23,87,90,91,92,97 };
    //int Tower[2] = { 88,89 };
    //for (int i = 33; i < 51; i++)
    //{
    //    landscape_wall[i] = 1;
    //}
    //for (int i = 0; i < NPoint; i++)
    //{
    //    for (int j = 0; j < 2; j++)
    //    {
    //        if (i == Tower[j])
    //            landscape_Tower[i] = 1;
    //    }
    //}
    //read0(arr, NPoint, MM, NN);
    //int** temploc = inquire_location0(arr, NPoint, MM, NN, landscape_wall);
    //file = fopen("F:/compute_verification/location/wall18.txt", "w"); // Open file for writing in text mode
    //if (file == NULL) {
    //    printf("Failed to open file!");
    //    return -1;
    //}
    //for (int i = 0; i < MM; i++)
    //{
    //    for (int j = 0; j < NN; j++)
    //        fprintf(file, "%u ", temploc[i][j]);
    //    fprintf(file, "\n");
    //}
    //temploc = inquire_location0(arr, NPoint, MM, NN, landscape_Tower);
    //file = fopen("F:/compute_verification/location/tower.txt", "w"); // Open file for writing in text mode
    //if (file == NULL) {
    //    printf("Failed to open file!");
    //    return -1;
    //}
    //for (int i = 0; i < MM; i++)
    //{
    //    for (int j = 0; j < NN; j++)
    //        fprintf(file, "%u ", temploc[i][j]);
    //    fprintf(file, "\n");
    //}
    //fclose(file);
    //for (int q = 0; q < MM; q++)
    //    delete (temploc[q]);
    //delete temploc;
    return 0;
}

/**
 * Compares two arrays for equality
 * @param a First array pointer
 * @param b Second array pointer
 * @param len Length of arrays to compare
 * @return true if arrays are identical, false otherwise
 */
bool compare(unsigned int* a, unsigned int* b, int len)
{
    for (int i = 0; i < len; i++)
        if (a[i] != b[i])
            return false;
    return true;
}

/**
 * Searches for a matching pattern in the Area structure
 * @param ar Pointer to Area structure containing patterns
 * @param data Pattern to search for
 * @return Pointer to matching Node if found, NULL otherwise
 */
Node* lookfor(Area* ar, unsigned int* data)
{
    Node* head = ar->point;
    while (head != NULL)
    {
        if (compare(head->arr, data, ar->len))
        {
            return head;
        }
        head = head->down;
    }
    return NULL;
}

/**
 * First compression stage: Bit-packing compression
 * Compresses binary landscape data by packing multiple layers into integers
 * @param x Original 3D binary array
 * @param len Number of layers
 * @param M Number of rows
 * @param N Number of columns
 * @return Compressed 3D array
 */
unsigned int*** compress1(unsigned int*** x, int len, int M, int N)
{
    int lenComp = len / 32 + (len % 32 == 0 ? 0 : 1);
    unsigned int*** p = (unsigned int***)malloc(lenComp * sizeof(unsigned int**));
    for (int i = 0; i < lenComp; i++)
    {
        p[i] = (unsigned int**)malloc(M * sizeof(unsigned int*));
        for (int j = 0; j < M; j++)
            p[i][j] = (unsigned int*)malloc(N * sizeof(unsigned int));
    }
    //for(int k=0;k<lenComp;k++){
    for (int row = 0; row < M; row++)
        for (int col = 0; col < N; col++)
        {
            unsigned int t = 0;
            int k = 0;
            for (int i = 0; i < lenComp; i++)//&& k < len
            {
                t = 0;
                for (int j = 0; j < 32; j++)
                {
                    //if (x[k][row][col] == true)
                    if (x[k][row][col] == true)
                    {
                        int js = 1;
                        for (int r = 1; r <= k % 32; r++)
                        {
                            js *= 2;
                        }
                        t += js;
                    }
                    k++;
                    if (k == len)
                        break;
                }
                p[i][row][col] = t;
            }
        }
    return p;
}

void write1(unsigned int*** goal, int len, int M, int N)
{
    FILE* file = fopen("f:/outVisibilityData/firstCompress.txt", "w");

    if (file == NULL) {
        printf("Failed to open file!");
        return;
    }
    for (int k = 0; k < len; k++) {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                fprintf(file, "%u ", goal[k][i][j]); // Write each element to the file, separated by spaces
            }
            fprintf(file, "\n"); // Write a newline character to indicate end of line
        }
        fprintf(file, "\n\n\n");
    }
    fclose(file);
}

Area* unique(unsigned int*** data, int len, int M, int N)
{
    unsigned int* temp = (unsigned int*)malloc(len * sizeof(unsigned int));
    Area* ha = (Area*)malloc(sizeof(Area));
    ha->len = len;
    ha->point = NULL;
    int size = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < len; k++)
                temp[k] = data[k][i][j];
            Node* p = lookfor(ha, temp);
            if (NULL == p)//Create a new landscape visual state
            {
                size++;
                Node* q = (Node*)malloc(sizeof(Node));
                q->arr = (unsigned int*)malloc(len * sizeof(unsigned int));
                for (int t = 0; t < len; t++)
                    q->arr[t] = temp[t];
                q->right = (coord*)malloc(sizeof(coord));
                q->right->r = i;
                q->right->c = j;
                q->right->next = NULL;
                q->down = ha->point;
                ha->point = q;
            }
            else //Add coordinates to the end of existing landscape visual state
            {
                coord* q = (coord*)malloc(sizeof(coord));
                q->r = i;
                q->c = j;
                q->next = p->right;
                p->right = q;
            }
        }
    }
    ha->size = size;
    return ha;
}

unsigned int** compress2(Area* comArea, unsigned int** data)
{
    unsigned int** data_status = (unsigned int**)malloc(comArea->size * sizeof(unsigned int*));
    for (int i = 0; i < comArea->size; i++)
    {
        data_status[i] = (unsigned int*)malloc(comArea->len * sizeof(unsigned int));
    }
    Node* l = comArea->point;
    int t = 0;
    while (l != NULL)
    {
        for (int k = 0; k < comArea->len; k++)
        {
            data_status[t][k] = l->arr[k];
        }
        coord* cp = l->right;
        while (cp != NULL)
        {
            data[cp->r][cp->c] = t;
            cp = cp->next;
        }
        l = l->down;
        t++;
    }
    return data_status;
}

void write2(Area* comArea, unsigned int** data_status, unsigned int** data, int M, int N)
{
    FILE* file = fopen("f:/outVisibilityData/secondCompress.txt", "w"); // Open file for writing in text mode
    if (file == NULL) {
        printf("Failed to open file!");
        return;
    }
    //Print landscape states
    for (int i = 0; i < comArea->size; i++)
    {
        fprintf(file, "%u:", i);
        for (int k = 0; k < comArea->len; k++)
            fprintf(file, "%u ", data_status[i][k]);
        fprintf(file, "\n");
    }
    //Print secondary compression data
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
            fprintf(file, "%u ", data[i][j]);
        fprintf(file, "\n");
    }
    fclose(file);
}

void write3(Area* comArea, unsigned int** data_status, unsigned int** a, int M, int N)
{
    FILE* file = fopen("f:/outVisibilityData/thirdCompress.txt", "w"); // Open file for writing in text mode

    if (file == NULL) {
        printf("Failed to open file!");
        return;
    }
    //Print landscape states
    for (int i = 0; i < comArea->size; i++)
    {
        fprintf(file, "%u:", i);
        for (int k = 0; k < comArea->len; k++)
            fprintf(file, "%u ", data_status[i][k]);
        fprintf(file, "\n");
    }
    //Print perception details (third compression of perception grid)
    for (int i = 0; i < M; i++)
    {
        unsigned int flag = a[i][0];
        int count = 0;
        for (int j = 1; j < N; j++)//Run-length encoding for each row
        {
            if (a[i][j] == flag)
                count++;
            else
            {
                fprintf(file, "%u:%u,", flag, count);
                count = 1;
                flag = a[i][j];
            }
        }
        fprintf(file, "%u:%u,", flag, count);
        fprintf(file, "\n");
    }
}
void read0(unsigned int***& goal, int len, int M, int N) {
    FILE* file = fopen("f:/outVisibilityData/benchmark.txt", "r");
    if (file == NULL) {
        printf("Failed to open file!");
        return;
    }

    for (int k = 0; k < len; k++) {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                fscanf(file, "%u", &goal[k][i][j]);
            }
        }
    }

    fclose(file);
}

void read1(unsigned int***& goal, int len, int M, int N) {
    FILE* file = fopen("f:/outVisibilityData/firstCompress.txt", "r");
    if (file == NULL) {
        printf("Failed to open file!");
        return;
    }

    for (int k = 0; k < len; k++) {
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                fscanf(file, "%u", &goal[k][i][j]);
            }
        }
    }

    fclose(file);
}

void read2(Area*& comArea, unsigned int**& data_status, unsigned int**& data, int M, int N) {
    FILE* file = fopen("f:/outVisibilityData/secondCompress.txt", "r");
    if (file == NULL) {
        printf("Failed to open file!");
        return;
    }

    // Read landscape states
    for (int i = 0; i < comArea->size; i++) {
        unsigned int index;
        fscanf(file, "%u:", &index);
        for (int k = 0; k < comArea->len; k++) {
            fscanf(file, "%u", &data_status[i][k]);
        }
    }

    // Read secondary compression data
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            fscanf(file, "%u", &data[i][j]);
        }
    }

    fclose(file);
}

void read3(Area*& comArea, unsigned int**& data_status, unsigned int**& a, int M, int N) {
    FILE* file = fopen("f:/outVisibilityData/thirdCompress.txt", "r");
    if (file == NULL) {
        printf("Failed to open file!");
        return;
    }

    // Read landscape states
    for (int i = 0; i < comArea->size; i++) {
        unsigned int index;
        fscanf(file, "%u:", &index);
        for (int k = 0; k < comArea->len; k++) {
            fscanf(file, "%u", &data_status[i][k]);
        }
    }

    // Read run-length encoding data
    for (int i = 0; i < M; i++) {
        unsigned int flag, count;
        int j = 0;
        while (fscanf(file, "%u:%u,", &flag, &count) == 2) {
            for (int k = 0; k < count + 1; k++) {
                if (j < N) {
                    a[i][j++] = flag;
                }
            }
        }
    }

    fclose(file);
}