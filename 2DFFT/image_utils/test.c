#include <stdio.h>
#include "1080p_soccer.h"
#include "img_write.h"
//#include "bilinear.h"
// #include <math.h>

#define WIDTH (unsigned int)1920
#define HEIGTH (unsigned int)1080
int main (void){

    img_write(imgIn, WIDTH, HEIGTH);
    return 0;
    
}
