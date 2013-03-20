#include <stdio.h>

void img_write_dat(unsigned int *imgData, unsigned int width, unsigned int height){

    unsigned int i = 0;
    unsigned int j = 0;

    FILE *fp;

    fp = fopen("img_in_1080p_16.dat", "w");

    /* fprintf(fp, "P2\n"); */
    /* fprintf(fp, "%d %d\n", width, height); */
    /* fprintf(fp, "255"); */

    /* for (i = 0; i < height; i++){ */
    /*     /\* fprintf(fp, "\n"); *\/ */
    /*     for (j = 0; j < width; j++){ */
    /*         fprintf(fp,"%d\n", imgData[i * width + j]); */
    /*     } */
    /* } */

    for (i = 0; i < height; i++){
        /* fprintf(fp, "\n"); */
        for (j = 0; j < width/2; j++){
            fprintf(fp,"%d\n", (imgData[i * width + 2*j]) | (imgData[i * width + 2*j+1]<<16));
        }
    }

    fclose(fp);
}

