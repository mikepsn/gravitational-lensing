#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN_GREY_VAL 0
#define MAX_GREY_VAL 255

#define G   6.67E-11    /* Gravitational Constant   */
#define C2  9.0E16      /* Speed of light squared   */

typedef struct
{
    int width;
    int height;
    int maxGreyVal;
    int **pixels;
} PGM;

PGM* readPGM(char* filename);
void outputPGM(PGM *pgm);
PGM* gravitationalLense(PGM* pgm, float D, float M);
void* memCalloc(size_t n, size_t size);
void* memAlloc(size_t size);

int main(int argc, char **argv)
{
    PGM* image = (PGM*)memAlloc(sizeof(PGM));
    float D = 0.0, M = 0.0;

    if (argc == 4)
    {
        image = readPGM(argv[1]);
        D = atof(argv[2]);
        M = atof(argv[3]);
        image = gravitationalLense(image, D, M);
        outputPGM(image);
    }
    else
    {
        fprintf(stderr,"Usage: image file.pgm D M\n");
    }
    return(0);
}

PGM* readPGM(char* filename)
{
    PGM* image = (PGM*)memAlloc(sizeof(PGM));
    char header[10];
    int i = 0, j = 0;
    FILE* fp;

    if ( (fp = fopen(filename, "r")) != (FILE*)NULL)
    {
        fscanf(fp, "%s", header);
        fscanf(fp, "%d %d", &(image->width), &(image->height));
        fscanf(fp, "%d", &(image->maxGreyVal));

        image->pixels = (int**)memCalloc(image->height, sizeof(int*));

        for (i = 0; i < image->height; i++)
        {
            image->pixels[i] = (int*)memCalloc(image->width, sizeof(int));

            for (j = 0; j < image->width; j++)
                fscanf(fp, "%d", &(image->pixels[i][j]));
        }

        fclose(fp);
    }
    else
    {
        fprintf(stderr,"Unable to open file\n");
        exit(1);
    }

    return(image);
}

PGM* gravitationalLense(PGM* pgm, float D, float M)
{
    PGM* image = (PGM*)memAlloc(sizeof(PGM));
    int i = 0, j = 0;
    float xc = 0.0, yc = 0.0;
    float xnew = 0.0, ynew = 0.0;
    float xp = 0.0, yp = 0.0;
    float xm = 0.0, ym = 0.0;
    int xi = 0, yi = 0;
    float beta = 0.0, psi = 0.0, phi = 0.0, angle = 0.0;
    float thetaE = 0.0, thetaPlus = 0.0, thetaMinus = 0.0;

    image->width = pgm->width;
    image->height = pgm->height;
    image->maxGreyVal = pgm->maxGreyVal;

    image->pixels = (int**)memCalloc(image->height, sizeof(int*));

    for (i = 0; i < image->height; i++)
    {
        image->pixels[i] = (int*)memCalloc(image->width, sizeof(int));

        for (j = 0; j < image->width; j++)
            image->pixels[i][j] = MIN_GREY_VAL;
    }
        /* Calculate center of image */

    xc = (image->width)/2.0; 
    yc = (image->height)/2.0;

        /* Calculate Einstein Ring Radius */

    thetaE = sqrt(4*G*M*D/C2);

        /* Calculate image position for each 
           pixel in the source */

    for (i = 0; i < image->height; i++)
    {
        for (j = 0; j < image->width; j++)
        {
            xnew =  j - xc ; /*j - xc;*/
            ynew = -i + yc;  /*-i + yc;*/

            beta = sqrt(xnew*xnew + ynew*ynew);

            if (xnew)
            {
                angle = atan(ynew/xnew);
            }
            else
            {
                if (ynew > 0.0)
                    angle = M_PI/2.0;
                else
                    angle = -M_PI/2.0;
            }
            angle = atan(ynew/xnew);

                /* Get correct angle for all quadrants */

            if (xnew > 0.0)
            {
                phi = angle;
                psi = angle + M_PI;
            }
            else if (xnew < 0.0)
            {
                phi = angle + M_PI;
                psi = angle;
            }
            
            if (xnew == 0.0)
            {
                if (ynew > 0.0)
                {
                    phi = M_PI/2.0;
                    psi = 3*M_PI/2.0;
                }
                else
                {
                    phi = 3*M_PI/2.0;
                    psi = M_PI/2.0;
                }
            }

            if (ynew == 0.0)
            {
                if (xnew > 0.0)
                {
                    phi = 0.0;
                    psi = M_PI;
                }
                else
                {
                    phi = M_PI;
                    psi = 0.0;
                }
            }

                /* Calculate thetaPlus and thetaMinus */

            thetaPlus  = 0.5*(beta + sqrt(beta*beta + 4*thetaE*thetaE));
            xp = thetaPlus*cos(phi);
            yp = thetaPlus*sin(phi);

            xi = (int)(xp + xc); 
            yi = -(int)(yp - yc); 
 
            if ((xi >= 0 && xi < image->width) && (yi >= 0 && yi < image->height))
            { 
                image->pixels[yi][xi] = pgm->pixels[i][j];
            }

            thetaMinus = 0.5*(beta - sqrt(beta*beta + 4*thetaE*thetaE));
            xm = thetaMinus*cos(psi);
            ym = thetaMinus*sin(psi);

            xi = (int)(xm + xc); 
            yi = -(int)(ym - yc); 

            if ((xi >= 0 && xi < image->width) && (yi >= 0 && yi < image->height))
            { 
                image->pixels[yi][xi] = pgm->pixels[i][j];
            }

        }
    }

    return(image);
}

void outputPGM(PGM *pgm)
{
    int i = 0, j = 0;
    int pixelCount = 0;

    printf("P2\n");
    printf("%d %d\n", pgm->width, pgm->height);
    printf("%d \n", pgm->maxGreyVal);

    for (i = 0; i < pgm->height; i++)
    {
        for (j = 0; j < pgm->width; j++)
        {
            printf("%d ", pgm->pixels[i][j]);
            if (pixelCount++ >= 10)
            {
                printf("\n");
                pixelCount = 0;
            }
        }

        if (pixelCount != 0)
        {
            printf("\n"); 
            pixelCount = 0;
        }
    }
}

void* memCalloc(size_t n, size_t size)
{
    void* tmp;

    tmp = calloc(n, size);

    if (tmp == (void*)NULL)
    {
        fprintf(stderr,"Unable to allocate memory...Exiting\n");
        exit(1);
    }

    return(tmp);
}

void* memAlloc(size_t size)
{
    return(memCalloc(1, size));
}


