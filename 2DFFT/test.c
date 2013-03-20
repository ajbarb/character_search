#include <stdio.h>
// #include "man.h"
// #include "img_write.h"

 #include <xdc/runtime/Timestamp.h>
#include "edma.h"

#define NUM_COLS 1024
#define NUM_ROWS 1024

/* #define L2SRAM_ADDR_SEG0 0x82000000 */
/* #define L2SRAM_ADDR_SEG1 0x83000000 */

#define NUM_BYTES 1024*32 -1
#define NUM_ARRAY 4
#define NUM_FRAMES 1
#define SRCBIDX NUM_BYTES
#define DSTBIDX NUM_BYTES

#define DMA_TX_SIZE (NUM_BYTES) * (NUM_ARRAY)
#define L2SRAM_ADDR_SEG_I0 0x00800000
#define L2SRAM_ADDR_SEG_I1 (L2SRAM_ADDR_SEG_I0) + (DMA_TX_SIZE)
#define L2SRAM_ADDR_SEG_O0 (L2SRAM_ADDR_SEG_I1) + (DMA_TX_SIZE)
#define L2SRAM_ADDR_SEG_O1 (L2SRAM_ADDR_SEG_O0) + (DMA_TX_SIZE)

// #define L2SRAM_ADDR_SEG1 0x00830000

/* #define NUM_BYTES 8 */
/* #define NUM_ARRAY 1024*32 - 1 */
/* #define NUM_FRAMES 1 */
/* #define SRCBIDX NUM_BYTES */
/* #define DSTBIDX NUM_BYTES */

int test (void){

    unsigned int i =0;
    unsigned int j =0;
    
    // const unsigned char *in;
    const unsigned char *in;
    unsigned char *out;

    unsigned int delta = 0;
    unsigned int timestamp = 0;

    unsigned char *l2Mem0;
    unsigned char *l2Mem1;

    unsigned char *outL2Mem0;
    unsigned char *outL2Mem1;

    unsigned char *ddr3Mem0;
    unsigned char *ddr3Mem1;

    in  = (unsigned char *)(0x80000000);
    out = (unsigned char *)(0x81000000);

    unsigned short  theta = 280;

    l2Mem0 = (unsigned char *)(L2SRAM_ADDR_SEG_I0);
    l2Mem1 = (unsigned char *)(L2SRAM_ADDR_SEG_I1);

    outL2Mem0 = (unsigned char *)(L2SRAM_ADDR_SEG_O0);
    outL2Mem1 = (unsigned char *)(L2SRAM_ADDR_SEG_O1);

    ddr3Mem0 = (unsigned char *)(0x84000000);
    ddr3Mem1 = (unsigned char *)(0x85000000);

    /* // generate horizontal gradient */
    /* for (i = 0; i < NUM_ROWS; i++) { */
    /*     // if (i > 511) pix = 255; */
    /*     // else pix = 0; */
    /*     for (j = 0; j < NUM_COLS; j++){ */
    /*         in[i*NUM_COLS + j] = (j) & (unsigned int)(0xFF); */
    /*     } */
    /* } */

    //check CPU copy from DDR3 to L2SRAM 
    /* for (i = 0; i < (NUM_ROWS)/2; i++) { */
    /*     // if (i > 511) pix = 255; */
    /*     // else pix = 0; */
    /*     for (j = 0; j < NUM_COLS; j++){ */
    /*         l2Mem0[i*NUM_COLS + j] = in[i*NUM_COLS + j]; */
    /*     } */
    /* } */
    
    edmaInit();
    timestamp = Timestamp_get32();

    edmaInitiateXfer( l2Mem0, in+0*(DMA_TX_SIZE), (int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES),
          SRCBIDX, DSTBIDX, 1,
          1, 0, 1);
    edmaWait4Completion(0);

    for (i = 0; i < 8; i++) {

        if ((i & (unsigned int)(0x01)) == 0){
            edmaInitiateXfer( l2Mem1, in+((i+1)*(DMA_TX_SIZE)), (int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES),
                              SRCBIDX, DSTBIDX, 1,
                              1, 0, 1);
            if (i != 0)
                edmaInitiateXfer(out+((i-1)*(DMA_TX_SIZE)), outL2Mem1,(int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES),
                  SRCBIDX, DSTBIDX, 1,
                  1, 1, 1);
            img_invert( l2Mem0, (NUM_COLS), (NUM_ROWS)/8, outL2Mem0);
            if (i != 0)
                 edmaWait4Completion(1);
            edmaWait4Completion(0);

        }
        else {
            edmaInitiateXfer( l2Mem0, in+((i+1)*(DMA_TX_SIZE)), (int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES),
                              SRCBIDX, DSTBIDX, 1,
                              1, 0, 1);
            edmaInitiateXfer(out+((i-1)*(DMA_TX_SIZE)), outL2Mem0,(int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES),
                             SRCBIDX, DSTBIDX, 1,
                             1, 1, 1);
            img_invert( l2Mem1, (NUM_COLS), (NUM_ROWS)/8, outL2Mem1);
            edmaWait4Completion(0);
            edmaWait4Completion(1);
        }
    }
    edmaInitiateXfer(out+((i-1)*(DMA_TX_SIZE)), outL2Mem1,(int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES),
                  SRCBIDX, DSTBIDX, 1,
                  1, 1, 1);
    edmaWait4Completion(1);
    /* edmaInitiateXfer( l2Mem1, in+3*(DMA_TX_SIZE), (int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES), */
    /*       SRCBIDX, DSTBIDX, 1, */
    /*       1, 0, 1); */
    /* edmaInitiateXfer(out+1*(DMA_TX_SIZE), outL2Mem1,(int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES), */
    /*       SRCBIDX, DSTBIDX, 1, */
    /*       1, 1, 1); */
    /* img_invert( l2Mem0, (NUM_COLS), (NUM_ROWS)/8, outL2Mem0); */
    /* edmaWait4Completion(0); */
    /* edmaWait4Completion(1); */

    /* edmaInitiateXfer( l2Mem0, in+4*(DMA_TX_SIZE), (int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES), */
    /*       SRCBIDX, DSTBIDX, 1, */
    /*       1, 0, 1); */
    /* edmaInitiateXfer(out+2*(DMA_TX_SIZE), outL2Mem0,(int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES), */
    /*       SRCBIDX, DSTBIDX, 1, */
    /*       1, 1, 1); */
    /* img_invert( l2Mem1, (NUM_COLS), (NUM_ROWS)/8, outL2Mem1); */
    /* edmaWait4Completion(0); */
    /* edmaWait4Completion(1); */

    /* edmaInitiateXfer(out+3*(DMA_TX_SIZE), outL2Mem1,(int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES), */
    /*       SRCBIDX, DSTBIDX, 1, */
    /*       1, 1, 1); */
    /* edmaWait4Completion(1); */
    

    /* edmaInitiateXfer( l2Mem1, in+(NUM_BYTES)*(NUM_ARRAY), (int)(NUM_BYTES), (int)(NUM_ARRAY), (int)(NUM_FRAMES), */
    /*       SRCBIDX, DSTBIDX, 1, */
    /*       1, 1, 1); */
    /* edmaWait4Completion(1); */

    edmaClose();

    // bilinear(in, (unsigned short)(NUM_COLS), (unsigned short) (NUM_ROWS), theta, out);
    /* unsigned int error =0; */
    /* // check transfer */
    /* for (i = 0; i < NUM_ROWS/4; i++) { */
    /*     for (j = 0; j < NUM_COLS; j++){ */
    /*         if (in[i*NUM_COLS + j] != l2Mem0[i*NUM_COLS + j]){ */
    /*             error += 1; */
    /*             /\* printf("in[%d]:%d != l2Mem0[%d]:%d\n", i*NUM_COLS +j, in[ i*NUM_COLS +j],  i*NUM_COLS +j, l2Mem0[ i*NUM_COLS +j]); *\/ */
    /*         } */
    /*     } */
    /* } */
    delta = Timestamp_get32() - timestamp;
    // printf("# of Clock Cycles: %d", *timestamp);
    printf("# of Clock Cycles: %d\n", delta);
    printf("Time: %fms\n", ((float)delta)/1000000);
    /* printf("Num DMA error = %d\n", error); */
    /* // img_write(out, 1024, 1024); */
    return 0;
    
}

