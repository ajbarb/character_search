#include <stdio.h>
#include"fft.h"

void img_write(mycomplex **imgData, int width, int height){

    unsigned int i = 0;
    unsigned int j = 0;

    FILE *fp;

    fp = fopen("img_out.pgm", "w");

    //fprintf(fp, "P2\n");
    //fprintf(fp, "%d %d\n", width, height);
    //fprintf(fp, "255");

    for (i = 0; i < height; i++){
        fprintf(fp, "\n");
        for (j = 0; j < width; j++){
            fprintf(fp, "%d ", (int)(imgData[i][j].r));
        }
    }

    fclose(fp);
}

