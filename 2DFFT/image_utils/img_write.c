#include <stdio.h>

void img_write(unsigned int *imgData, unsigned int width, unsigned int height){

    unsigned int i = 0;
    unsigned int j = 0;

    FILE *fp;

    fp = fopen("img_out.pgm", "w");

    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255");

    for (i = 0; i < height; i++){
        fprintf(fp, "\n");
        for (j = 0; j < width; j++){
            fprintf(fp, "%d ", imgData[i * width + j]);
        }
    }

    fclose(fp);
}

