//SuperPixel.c


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include "mex.h"
#include <semaphore.h>

#ifndef max
#define max(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef min
#define min(a, b) ((a)<(b)?(a):(b))
#endif





float dist_bresenham_contour(int x1,int y1,int xc,int yc, float * contour, int height)
{
    int dx, dy, i, e;
    int incx, incy, inc1, inc2;
    int x,y;
    float count = 0.0;
    float max_count = 20;
    float max_contour = 0;
    
    dx = xc - x1;
    dy = yc - y1;
    
    if(dx < 0)
        dx = -dx;
    if(dy < 0)
        dy = -dy;
    
    incx = 1;
    if(xc < x1)
        incx = -1;
    
    incy = 1;
    if(yc < y1)
        incy = -1;
    
    x=x1;
    y=y1;
    
    int ind;
    
    
//     count += 1;
    if(dx > dy)
    {
        //draw_pixel(x,y, BLACK);
        e = 2*dy - dx;
        inc1 = 2*( dy -dx);
        inc2 = 2*dy;
        for(i = 0; i < dx; i++)
        {
            if(e >= 0)
            {
                y += incy;
                e += inc1;
            }
            else
                e += inc2;
            
            x += incx;
            //count += 1;
            ind = x*height + y;
            
            if (count > max_count)
                return max_contour;
//                 return dist/count;
            else {
//                 if (klabels[ind] == n) {
//                 dist += contour[ind];
                if (contour[ind] > max_contour)
                    max_contour = contour[ind];
                count += 1;
//                 }
            }
            
            
            //dist += contour[ind];
            
        }
    }
    else
    {
        //draw_pixel(x,y, BLACK);
        e = 2*dx - dy;
        inc1 = 2*( dx - dy);
        inc2 = 2*dx;
        
        for(i = 0; i < dy; i++)
        {
            if(e >= 0)
            {
                x += incx;
                e += inc1;
            }
            else
                e += inc2;
            
            y += incy;
            //count += 1;
            
            ind = x*height + y;
            
            if (count > max_count)
                return max_contour;
//                 return dist/count;
            else {
//              if (klabels[ind] == n) {
//                 dist += contour[ind];
                if (contour[ind] > max_contour)
                    max_contour = contour[ind];
                count += 1;
//              }
            }
        }
    }
//     if (dist/count > 2)
//         printf("%lf\n", dist/count);
    
//     if (count > 0)
//         dist /= count;
    
//     return dist;
    
    return max_contour;
}






float dist_bresenham(int x1,int y1,int xc,int yc, float * m_lvec, float* m_avec, float *m_bvec,
        float lc, float ac, float bc, int height)
{
    int dx, dy, i, e;
    int incx, incy, inc1, inc2;
    int x,y;
    float count = 0.0;
    float max_count = 20;
    
    dx = xc - x1;
    dy = yc - y1;
    
    if(dx < 0)
        dx = -dx;
    if(dy < 0)
        dy = -dy;
    
    incx = 1;
    if(xc < x1)
        incx = -1;
    
    incy = 1;
    if(yc < y1)
        incy = -1;
    
    x=x1;
    y=y1;
    
    float dist = 0;
    int ind;
    
    
//     count += 1;
    if(dx > dy)
    {
        //draw_pixel(x,y, BLACK);
        e = 2*dy - dx;
        inc1 = 2*( dy -dx);
        inc2 = 2*dy;
        for(i = 0; i < dx; i++)
        {
            if(e >= 0)
            {
                y += incy;
                e += inc1;
            }
            else
                e += inc2;
            
            x += incx;
            //count += 1;
            
            ind = x*height + y;
            
            if (count > max_count)
                return dist/count;
            else {
                
//                 if (klabels[ind] == n) {
                dist +=  (lc-m_lvec[ind])*(lc-m_lvec[ind])
                +(ac-m_avec[ind])*(ac-m_avec[ind])
                +(bc-m_bvec[ind])*(bc-m_bvec[ind]);
                count += 1;
//                 }
            }
        }
    }
    else
    {
        //draw_pixel(x,y, BLACK);
        e = 2*dx - dy;
        inc1 = 2*( dx - dy);
        inc2 = 2*dx;
        
        for(i = 0; i < dy; i++)
        {
            if(e >= 0)
            {
                x += incx;
                e += inc1;
            }
            else
                e += inc2;
            
            y += incy;
            //count += 1;
            
            ind = x*height + y;
            if (count > max_count)
                return dist/count;
            else {
//                 if (klabels[ind] == n) {
                dist +=  (lc-m_lvec[ind])*(lc-m_lvec[ind])
                +(ac-m_avec[ind])*(ac-m_avec[ind])
                +(bc-m_bvec[ind])*(bc-m_bvec[ind]);
                count += 1;
//                 }
                
            }
        }
    }
    
//     printf("%lf\n", dist/count);
    if (count > 0)
        return dist/count;
    else
        return dist;
    
}





void RGB2XYZ(
        const float		sR,
        const float		sG,
        const float		sB,
        float*			X,
        float*			Y,
        float*			Z)
{
//     float R = sR/255.0;
//     float G = sG/255.0;
//     float B = sB/255.0;
    
    float R = sR;
    float G = sG;
    float B = sB;
    
    float r, g, b;
    
    if (R <= 0.04045)
        r = R/12.92;
    else
        r = pow((R+0.055)/1.055,2.4);
    
    if (G <= 0.04045)
        g = G/12.92;
    else
        g = pow((G+0.055)/1.055,2.4);
    
    if (B <= 0.04045)
        b = B/12.92;
    else
        b = pow((B+0.055)/1.055,2.4);
    
    
    *X = r*0.4124564 + g*0.3575761 + b*0.1804375;
    *Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
    *Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
    
}



void RGB2LAB(const float sR, const float sG, const float sB, float* lval, float* aval, float* bval)
{
    //------------------------
    // sRGB to XYZ conversion
    //------------------------
    float X, Y, Z;
    RGB2XYZ(sR, sG, sB, &X, &Y, &Z);
    
    //------------------------
    // XYZ to LAB conversion
    //------------------------
    float epsilon = 0.008856;	//actual CIE standard
    float kappa   = 903.3;		//actual CIE standard
    
    float Xr = 0.950456;	//reference white
    float Yr = 1.0;		//reference white
    float Zr = 1.088754;	//reference white
    
    float xr = X/Xr;
    float yr = Y/Yr;
    float zr = Z/Zr;
    
    float fx, fy, fz;
    if(xr > epsilon)
        fx = pow(xr, 1.0/3.0);
    else
        fx = (kappa*xr + 16.0)/116.0;
    
    if(yr > epsilon)
        fy = pow(yr, 1.0/3.0);
    else
        fy = (kappa*yr + 16.0)/116.0;
    
    if(zr > epsilon)
        fz = pow(zr, 1.0/3.0);
    else
        fz = (kappa*zr + 16.0)/116.0;
    
    *lval = 116.0*fy-16.0;
    *aval = 500.0*(fx-fy);
    *bval = 200.0*(fy-fz);
}


void DoRGBtoLABConversion(
        float*      ubuffr,
        float*       ubuffg,
        float*       ubuffb,
        float*      lvec,
        float*      avec,
        float*      bvec,
        const int   m_height,
        const int   m_width)
{
    
    int sz = m_width*m_height;
    
    for( int j = 0; j < sz; j++ )
    {
        float r = ubuffr[j];
        float g = ubuffg[j];
        float b = ubuffb[j];
        
        RGB2LAB(r, g, b, &lvec[j], &avec[j], &bvec[j]);
    }
    
}


void DetectLabEdges(
        const float*				lvec,
        const float*				avec,
        const float*				bvec,
        const int					height,
        const int					width,
        float*     				edges)
{
    for( int j = 1; j < height-1; j++ )
    {
        for( int k = 1; k < width-1; k++ )
        {
            int i = j+k*height;
            
            float dx = (lvec[i-1]-lvec[i+1])*(lvec[i-1]-lvec[i+1]) +
                    (avec[i-1]-avec[i+1])*(avec[i-1]-avec[i+1]) +
                    (bvec[i-1]-bvec[i+1])*(bvec[i-1]-bvec[i+1]);
            
            float dy = (lvec[i-width]-lvec[i+width])*(lvec[i-width]-lvec[i+width]) +
                    (avec[i-width]-avec[i+width])*(avec[i-width]-avec[i+width]) +
                    (bvec[i-width]-bvec[i+width])*(bvec[i-width]-bvec[i+width]);
            
            edges[i] = dx*dx + dy*dy;
        }
    }
}








void GetLABXYSeeds_ForGivenStepSize(
        float *               m_lvec,
        float *               m_avec,
        float *               m_bvec,
        float *        		kseedsl,
        float *				kseedsa,
        float *				kseedsb,
        float *				kseedsx,
        float *				kseedsy,
        const int				STEP,
        const int               m_height,
        const int               m_width,
        const int               xstrips,
        const int               ystrips,
        const float            xerrperstrip,
        const float            yerrperstrip,
        const int				perturbseeds)
{
    const int hexgrid = 0;
    int numseeds = 0;
    int n = 0;
    
    int xoff = STEP/2;
    int yoff = STEP/2;
    //-------------------------
    numseeds = xstrips*ystrips;
    //-------------------------
    
    for( int y = 0; y < ystrips; y++ )
    {
        int ye = y*yerrperstrip;
        for( int x = 0; x < xstrips; x++ )
        {
            int xe = x*xerrperstrip;
            int seedx = (x*STEP+xoff+xe);
            
            if (hexgrid){
                seedx = x*STEP+(xoff<<(y&0x1))+xe;
                seedx = min(m_width-1,seedx);
            }//for hex grid sampling
            
            int seedy = (y*STEP+yoff+ye);
            int i = seedx*m_height + seedy;
            
            kseedsl[n] = m_lvec[i];
            kseedsa[n] = m_avec[i];
            kseedsb[n] = m_bvec[i];
            kseedsx[n] = seedx;
            kseedsy[n] = seedy;
            n++;
        }
    }
    
    
    
    
}


typedef struct{
    float * m_lvec;
    float * m_avec;
    float * m_bvec;
    float * kseedsl;
    float * kseedsa;
    float * kseedsb;
    float * kseedsx;
    float * kseedsy;
    int * klabels;
    int m_height;
    int m_width;
    float * distvec;
    int offset;
    int ini;
    int fin;
    float invwt;
    sem_t* mutex;
    float * kseedsx_p;
    float * kseedsy_p;
    float * maxlab;
    float * distlab;
    int bres_lab;
    int bres_contour;
    float * contour;
    float * Ip2;
    int pw;
    float * Ip;
    float * sp_hist;
    int h_bins;
    int itr;
    int numk;
    float * m_lvec01;
    float * m_avec01;
    float * m_bvec01;
}slic_mt_struct;



void *SLIC_core(void *arg) {
    
    slic_mt_struct inputs;
    inputs = *(slic_mt_struct*) arg;
    
    float * m_lvec = inputs.m_lvec;
    float * m_avec = inputs.m_avec;
    float * m_bvec = inputs.m_bvec;
    float * kseedsl  = inputs.kseedsl;
    float * kseedsa  = inputs.kseedsa;
    float * kseedsb  = inputs.kseedsb;
    float * kseedsx  = inputs.kseedsx;
    float * kseedsy  = inputs.kseedsy;
    int * klabels  = inputs.klabels;
    int m_height = inputs.m_height;
    int m_width = inputs.m_width;
    float * distvec = inputs.distvec;
    int offset = inputs.offset;
    float invwt = inputs.invwt;
    sem_t* mutex = inputs.mutex;
    float * kseedsx_p = inputs.kseedsx_p;
    float * kseedsy_p = inputs.kseedsy_p;
    int bres_lab = inputs.bres_lab;
    int bres_contour = inputs.bres_contour;
    float * contour = inputs.contour;
    float * Ip2 = inputs.Ip2;
    float * Ip = inputs.Ip;
    int pw = inputs.pw;
    float * sp_hist = inputs.sp_hist;
    int h_bins = inputs.h_bins;
    int itr = inputs.itr;
    int numk = inputs.numk;
    float * m_lvec01 = inputs.m_lvec01;
    float * m_avec01 = inputs.m_avec01;
    float * m_bvec01 = inputs.m_bvec01;
    
    int ini = inputs.ini;
    int fin = inputs.fin;
    
    int x1, y1, x2, y2;
    float l, a, b;
    float dist;
    float distxy, tmp;
    float lambda = 1;
    if (bres_lab > 0)
        lambda = 0.5;
    float avg_contour;
    float *avg_ptr = &avg_contour;
    *avg_ptr = 0.0;
    float dist_tmp;
    int sz = m_height*m_width;
    
    if ((itr == 0) && (bres_lab == 2))
        invwt *= 1000;
    
    for( int n = ini; n < fin; n++ )
    {
        y1 = max(0.0,			kseedsy[n]-offset);
        y2 = min((float)m_height,	kseedsy[n]+offset);
        x1 = max(0.0,			kseedsx[n]-offset);
        x2 = min((float)m_width,	kseedsx[n]+offset);
        
        
        for( int y = y1; y < y2; y++ )
        {
            for( int x = x1; x < x2; x++ )
            {
                int i = x*m_height + y;
                
//                 l = m_lvec[i];
//                 a = m_avec[i];
//                 b = m_bvec[i];
//                 dist =	lambda*((l - kseedsl[n])*(l - kseedsl[n]) +
//                         (a - kseedsa[n])*(a - kseedsa[n]) +
//                         (b - kseedsb[n])*(b - kseedsb[n]));
                
                
//                 dist = 0;
//                 float count = 0;
//                 for (int dx=-pw; dx<=pw; dx++) {
//                     for (int dy=-pw; dy<=pw; dy++) {
//                         if ((x+dx<m_width)&&(x+dx>=0)&(y+dy>=0)&(y+dy<m_height)){
//                             int pos = y+dy + (x+dx)*m_height;
//                             count += 1;
//                             dist += (m_lvec[pos] - kseedsl[n])*(m_lvec[pos] - kseedsl[n]) +
//                                     (m_avec[pos] - kseedsa[n])*(m_avec[pos] - kseedsa[n]) +
//                                     (m_bvec[pos] - kseedsb[n])*(m_bvec[pos] - kseedsb[n]);
//                         }
//                     }
//                 }
//                 dist = lambda*dist/count;
//
                //Speed up patch-based distance
                dist = Ip2[i] + kseedsl[n]*kseedsl[n] - 2*kseedsl[n]*Ip[i];
                dist += Ip2[i+sz] + kseedsa[n]*kseedsa[n] - 2*kseedsa[n]*Ip[i+sz];
                dist += Ip2[i+sz*2] + kseedsb[n]*kseedsb[n] - 2*kseedsb[n]*Ip[i+sz*2];
                dist *= lambda;
                
                
                if (bres_lab == 1) {
//                     if (sqrt( (x-(int) round(kseedsx_p[n]))*(x-(int) round(kseedsx_p[n])) + (y-(int) round(kseedsy_p[n]))*(y-(int) round(kseedsy_p[n])) ) <  5)
//                        dist += (1-lambda)*dist_bresenham(x,y,(int) round(kseedsx_p[n]),(int) round(kseedsy_p[n]),m_lvec,m_avec,m_bvec,kseedsl[n],kseedsa[n],kseedsb[n],m_height);
//                     else
                        dist += (1-lambda)*dist_bresenham(x,y,(int) round(kseedsx_p[n]),(int) round(kseedsy_p[n]),m_lvec,m_avec,m_bvec,kseedsl[n],kseedsa[n],kseedsb[n],m_height);
                }
//                 if ((bres_lab == 2) && (itr>0)) {
//                     float ratio_hist = 100;
//                     float dist_tmp = ratio_hist*dist_bresenham_histogram(x,y,(int) round(kseedsx_p[n]),(int) round(kseedsy_p[n]),m_lvec01,m_avec01,m_bvec01,sp_hist,h_bins,n,numk,m_height);
// //                     printf("%f, %f\n", dist, dist_tmp);
// //                     mexEvalString("drawnow");
//                     dist += dist_tmp;
//                 }
                
                if (bres_contour) {
//                     printf("hello\n");
//                     mexEvalString("drawnow");
                    tmp = dist_bresenham_contour(x,y,(int) round(kseedsx_p[n]), (int) round(kseedsy_p[n]),contour,m_height);
//                     printf("%lf\n", tmp);
//                     mexEvalString("drawnow");
//                     dist *= tmp;
//                     dist *= 1;
                    //avg_contour = 0.0;
//                     dist_tmp = (1-lambda)*dist_bresenham_combined(x,y,(int) round(kseedsx_p[n]),(int) round(kseedsy_p[n]),m_lvec,m_avec,m_bvec,
//                             kseedsl[n],kseedsa[n],kseedsb[n],m_height, n, klabels,contour); //, avg_ptr);
//                     dist = dist + dist_tmp;
                    
                }
                else {
                    tmp = 1;
                }
//                 printf("hello2\n");
//                     mexEvalString("drawnow");
                
                distxy =	(x - kseedsx[n])*(x - kseedsx[n]) +
                            (y - kseedsy[n])*(y - kseedsy[n]);

                //------------------------------------------------------------------------
                dist += distxy*invwt;//dist = sqrt(dist) + sqrt(distxy*invwt);//this is more exact
                //------------------------------------------------------------------------
                
                dist = dist*tmp;
                
                
                sem_wait(&(mutex[i]));
                if( dist < distvec[i] )  {
                    distvec[i] = dist;
                    klabels[i]  = n;
                }
                sem_post(&(mutex[i]));
                
            }
        }
    }
    return NULL;
    
}



void SLIC(
        float * m_lvec,
        float * m_avec,
        float * m_bvec,
        float *				kseedsl,
        float *				kseedsa,
        float *				kseedsb,
        float *				kseedsx,
        float *				kseedsy,
        int *					klabels,
        const int				STEP,
        const float				M,
        const int numseeds,
        const int m_height,
        const int m_width,
        int nb_iter,
        int bres_lab,
        int update_center,
        int bres_contour,
        float * contour,
        int pw,
        float * Ip2,
        float * Ip,
        int h_bins,
        float * m_lvec01,
        float * m_avec01,
        float * m_bvec01)
{
    int sz = m_width*m_height;
    const int numk = numseeds;
    //----------------
    int offset = STEP;
    
    float * clustersize = (float *)calloc(numk, sizeof(float));
    float * inv = (float *)calloc(numk, sizeof(float));//to store 1/clustersize[k] values
    
    float * sigmal = (float *)calloc(numk, sizeof(float));
    float * sigmaa = (float *)calloc(numk, sizeof(float));
    float * sigmab = (float *)calloc(numk, sizeof(float));
    float * sigmax = (float *)calloc(numk, sizeof(float));
    float * sigmay = (float *)calloc(numk, sizeof(float));
    float * kseedsx_p = (float *)calloc(numk, sizeof(float));
    float * kseedsy_p = (float *)calloc(numk, sizeof(float));
    
    float * min_dist_center = (float *)calloc(numk, sizeof(float));
    
    float * sp_hist = (float *)calloc(numk*h_bins*3, sizeof(float));
    float * distvec = (float *)calloc(sz, sizeof(float));
//     for (int i=0; i<sz; i++) {
//         distvec[i] = DBL_MAX;1
//     }
    
    for( int k = 0; k < numk; k++ )
    {
        kseedsx_p[k] = kseedsx[k];
        kseedsy_p[k] = kseedsy[k];
    }
    
    
    float invwt = 1.0/((STEP/M)*(STEP/M));
    
    
    int thread_nbr = 1; //get_nprocs_conf();
    if (thread_nbr == 0)
        thread_nbr = 1;
//     printf("thread nbr %d\n", thread_nbr);
    int step = floor(numk/thread_nbr);
    int * slice_vect = (int *) calloc(thread_nbr+1, sizeof(int));
    
    slice_vect[0] = 0;
    for (int i=1; i<thread_nbr; i++)
        slice_vect[i] = i*step;
    slice_vect[thread_nbr] = numk;
    
    
    //Mutex
    sem_t* mutex= (sem_t*) malloc(m_height*m_width*sizeof(sem_t));
    for (int i=0;i<m_height*m_width;i++)    {
        sem_init(&(mutex[i]), 0, 1);
    }
    
    //THREAD BUILDING
    pthread_t *thread_list_th = (pthread_t *) calloc(thread_nbr,
            sizeof(pthread_t));
    slic_mt_struct* ThreadArgs =
            (slic_mt_struct*)calloc(thread_nbr, sizeof(slic_mt_struct));
    
    for( int itr = 0; itr < nb_iter; itr++ )
    {
        
        for (int i=0; i<sz; i++) {
            distvec[i] = DBL_MAX;
        }
        
        /*Launching of the THREADS*/
        for (int i=0; i < thread_nbr; i++) {
            /*Thread arguments*/
            ThreadArgs[i].kseedsx = kseedsx;
            ThreadArgs[i].kseedsy = kseedsy;
            ThreadArgs[i].m_lvec = m_lvec;
            ThreadArgs[i].m_avec = m_avec;
            ThreadArgs[i].m_bvec = m_bvec;
            ThreadArgs[i].m_height = m_height;
            ThreadArgs[i].m_width = m_width;
            ThreadArgs[i].kseedsl = kseedsl;
            ThreadArgs[i].kseedsa = kseedsa;
            ThreadArgs[i].kseedsb = kseedsb;
            ThreadArgs[i].distvec = distvec;
            ThreadArgs[i].klabels = klabels;
            ThreadArgs[i].invwt = invwt;
            ThreadArgs[i].offset = offset;
            ThreadArgs[i].mutex = mutex;
            ThreadArgs[i].kseedsx_p = kseedsx_p;
            ThreadArgs[i].kseedsy_p = kseedsy_p;
            ThreadArgs[i].contour = contour;
            ThreadArgs[i].bres_lab = bres_lab;
            ThreadArgs[i].bres_contour = bres_contour;
            ThreadArgs[i].Ip = Ip;
            ThreadArgs[i].Ip2 = Ip2;
            ThreadArgs[i].pw = pw;
            ThreadArgs[i].sp_hist = sp_hist;
            ThreadArgs[i].h_bins = h_bins;
            ThreadArgs[i].numk = numk;
            ThreadArgs[i].itr = itr;
            ThreadArgs[i].m_lvec01 = m_lvec01;
            ThreadArgs[i].m_avec01 = m_avec01;
            ThreadArgs[i].m_bvec01 = m_bvec01;
            
            ThreadArgs[i].ini = slice_vect[i];
            ThreadArgs[i].fin = slice_vect[i+1];
            
            if (pthread_create(&thread_list_th[i], NULL,
                    SLIC_core, &ThreadArgs[i]))
                printf("Error creating a thread!\n");
        }
        
//         printf("exit\n");
//         mexEvalString("drawnow");
        
        
        
        /*Wait for all threads to end*/
        for (int j=0; j<thread_nbr; j++) {
            pthread_join(thread_list_th[j],NULL);
        }
        
        if (itr < nb_iter-1) {
            
            for(int i=0; i<numk; i++) {
                sigmal[i] = (float) 0;
                sigmaa[i] = (float) 0;
                sigmab[i] = (float) 0;
                sigmax[i] = (float) 0;
                sigmay[i] = (float) 0;
                clustersize[i] = (float) 0;
                if (bres_lab == 2) {
                    for (int l=0; l< h_bins*3; l++)
                        sp_hist[i+l*numk] = 0;
                }
            }
            
            int ind = 0;
            for( int r = 0; r < m_width; r++ )
            {
                for( int c = 0; c < m_height; c++ )
                {
                    int n = klabels[ind];
                    sigmal[n] += m_lvec[ind];
                    sigmaa[n] += m_avec[ind];
                    sigmab[n] += m_bvec[ind];
                    sigmax[n] += r;
                    sigmay[n] += c;
                    //------------------------------------
                    //edgesum[klabels[ind]] += edgemag[ind];
                    //------------------------------------
                    clustersize[n] += 1.0;
                    
                    if (bres_lab == 2) {
                        //hist filling
                        sp_hist[n + (int) min(floor(m_lvec01[ind]*h_bins),h_bins-1)*numk] += 1;
                        sp_hist[n + (int) min(floor(m_avec01[ind]*h_bins),h_bins-1)*numk + h_bins*numk] += 1;
                        sp_hist[n + (int) min(floor(m_bvec01[ind]*h_bins),h_bins-1)*numk + h_bins*numk*2] += 1;
                    }
                    
                    ind++;
                }
            }
            
            for( int k = 0; k < numk; k++ )
            {
                if( clustersize[k] <= 0 ) clustersize[k] = 1;
                inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
            }
            
            for( int k = 0; k < numk; k++ )
            {
                if(clustersize[k] > 0) {
                    kseedsl[k] = sigmal[k]*inv[k];
                    kseedsa[k] = sigmaa[k]*inv[k];
                    kseedsb[k] = sigmab[k]*inv[k];
                    kseedsx[k] = sigmax[k]*inv[k];
                    kseedsy[k] = sigmay[k]*inv[k];
                    //
                    kseedsx_p[k] = kseedsx[k];
                    kseedsy_p[k] = kseedsy[k];
                    
                    
                    //Cumulated histogram
                    if (bres_lab == 2) {
                        sp_hist[k] *= inv[k];
                        sp_hist[k+h_bins*numk] *= inv[k];
                        sp_hist[k+h_bins*2*numk] *= inv[k];
                        
                        for (int i=1; i<h_bins; i++)
                            sp_hist[k + i*numk] = sp_hist[k + i*numk]*inv[k] + sp_hist[k + (i-1)*numk] ;
                        
                        for (int i=h_bins+1; i<2*h_bins; i++)
                            sp_hist[k + i*numk] = sp_hist[k + i*numk]*inv[k] + sp_hist[k + (i-1)*numk] ;
                        
                        for (int i=2*h_bins+1; i<3*h_bins; i++)
                            sp_hist[k + i*numk] = sp_hist[k + i*numk]*inv[k] + sp_hist[k + (i-1)*numk] ;
                    }
                }
                
                else {
                    
                    kseedsx[k] = 9999;
                    kseedsy[k] = 9999;
                    
                    if (bres_lab == 2) {
                        for (int i=0; i<h_bins*3; i++)
                            sp_hist[k+i*numk] = 10000;
                    }
                }
                
            }
            
            if (update_center) {
                
                for( int k = 0; k < numk; k++ )  {
                    min_dist_center[k] = FLT_MAX;
                }
                
                //Update center position - projection to closest pixel
                ind = 0;
                float dist_center = 0;
                for( int r = 0; r < m_width; r++ )
                {
                    for( int c = 0; c < m_height; c++ )
                    {
                        int n = klabels[ind];
                        if (klabels[(int) round(kseedsy[n])+ (int) round(kseedsx[n])*m_height] != n) {
                            dist_center = (float) (kseedsx[n]-r)*(kseedsx[n]-r) + (kseedsx[n]-c)*(kseedsx[n]-c);
                            if (dist_center < min_dist_center[n]) {
                                min_dist_center[n] = dist_center;
                                kseedsy_p[n] = c;
                                kseedsx_p[n] = r;
                            }
                        }
                        
                        ind++;
                        
                    }
                }
                
            }
        }
        
        
    }
    
    free(sp_hist);
    free(clustersize);
    free(inv);
    free(distvec);
    free(sigmal);
    free(sigmaa);
    free(sigmab);
    free(sigmax);
    free(sigmay);
    free(min_dist_center);
    free(kseedsx_p);
    free(kseedsy_p);
    free(mutex);
    free(slice_vect);
    free(thread_list_th);
    free(ThreadArgs);
    
}



void EnforceLabelConnectivity(
        const int*					labels,//input labels that need to be corrected to remove stray labels
        const int					height,
        const int					width,
        int*						nlabels,//new labels
        int*						numlabels,//the number of labels changes in the end if segments are removed
        const int					K) //the number of superpixels desired by the user
{
//	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
//	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
    
    const int dx4[4] = {-1,  0,  1,  0};
    const int dy4[4] = { 0, -1,  0,  1};
    
    const int sz = width*height;
    const int SUPSZ = sz/(K);  //!!!!
    //nlabels.resize(sz, -1);
    for( int i = 0; i < sz; i++ )
        nlabels[i] = -1;
    
    int label = 0;
    
    int * xvec= (int *) calloc(sz, sizeof(int));
    int * yvec= (int *) calloc(sz, sizeof(int));
    
    int oindex = 0;
    int adjlabel = 0;//adjacent label
    for( int j = 0; j < width; j++ )
    {
        for( int k = 0; k < height; k++ )
        {
            if( 0 > nlabels[oindex] )
            {
                nlabels[oindex] = label;
                //--------------------
                // Start a new segment
                //--------------------
                xvec[0] = j;
                yvec[0] = k;
                //-------------------------------------------------------
                // Quickly find an adjacent label for use later if needed
                //-------------------------------------------------------
                for( int n = 0; n < 4; n++ )
                {
                    int x = xvec[0] + dx4[n];
                    int y = yvec[0] + dy4[n];
                    if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                    {
                        int nindex = x*height + y;
                        if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
                    }
                }
                
                int count = 1;
                for( int c = 0; c < count; c++ )
                {
                    for( int n = 0; n < 4; n++ )
                    {
                        int x = xvec[c] + dx4[n];
                        int y = yvec[c] + dy4[n];
                        
                        if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        {
                            int nindex = x*height + y;
                            
                            if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
                            {
                                xvec[count] = x;
                                yvec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }
                        
                    }
                }
                //-------------------------------------------------------
                // If segment size is less then a limit, assign an
                // adjacent label found before,
                //and decrement label count.
                //-------------------------------------------------------
                if(count <= SUPSZ >> 2)
                {
                    for( int c = 0; c < count; c++ )
                    {
                        int ind = yvec[c]+xvec[c]*height;
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }
    *numlabels = label;
    
//
    free(xvec);
    free(yvec);
}






// typedef struct{ for (int i=0; i< sz; i++) {
//     float * ubuffr;
//     float * ubuffg;
//     float * ubuffb;
//     int * klabels;
//
//     int h_img;
//     int w_img;
//     int offset;
//
//     int SP_nbr;
//     int superpixelsize;
//     float compactness;
//     int perturbseeds;
//
//     int nb_iter;
//     int bres_lab;
//     int update_center;
//     int bres_contour;
//     float * contour;
//
// }do_slic_mt_struct;


//
// void *DoSuperpixelSegmentation_ForGivenSuperpixelSize_thread(void *arg) {
//
//
//     do_slic_mt_struct inputs;
//     inputs = *(do_slic_mt_struct*) arg;
//
//     float * ubuffr = inputs.ubuffr;
//     float * ubuffg = inputs.ubuffg;
//     float * ubuffb = inputs.ubuffb;
//     int * klabels  = inputs.klabels;
//
//     int height = inputs.h_img;
//     int width = inputs.w_img;
//     int offset = inputs.offset;
//
//     int numlabels = inputs.SP_nbr;
//     int superpixelsize = inputs.superpixelsize;
//     const float compactness = inputs.compactness;
//     const int perturbseeds = inputs.perturbseeds;
//
//     int bres_lab =  inputs.bres_lab;
//     int update_center =  inputs.update_center;
//     int bres_contour = inputs.bres_contour;
//     float * contour = inputs.contour;
//
//     int nb_iter = inputs.nb_iter;
//
//
//
//     const int STEP = sqrt((float)superpixelsize)+0.5;
//     int m_width  = width;
//     int m_height = height;
//     int sz = m_width*m_height;
//
//     float* m_lvec = (float *)calloc(sz, sizeof(float));
//     float* m_avec = (float *)calloc(sz, sizeof(float));
//     float* m_bvec = (float *)calloc(sz, sizeof(float));
//
//     //--------------------------------------------------
//     DoRGBtoLABConversion(ubuffr, ubuffg, ubuffb, m_lvec, m_avec, m_bvec, m_height, m_width);
//     //--------------------------------------------------
//
//     float * edgemap = (float *)calloc(sz, sizeof(float));
//
//     //perturb seeds is not absolutely necessary, one can set this flag to false
//     if(perturbseeds)
//         DetectLabEdges(m_lvec, m_avec, m_bvec, m_height, m_width, edgemap);
//
//     int xstrips = 0.5+(float)m_width/(float)STEP;
//     int ystrips = 0.5+(float)m_height/(float)STEP;
//     int numseeds = xstrips*ystrips;
//
//     int xerr = m_width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = m_width - STEP*xstrips;}
//     int yerr = m_height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = m_height- STEP*ystrips;}
//
//     float xerrperstrip = (float)xerr/(float)xstrips;
//     float yerrperstrip = (float)yerr/(float)ystrips;
//
//     float * kseedsl = (float *)calloc(numseeds, sizeof(float));
//     float * kseedsa = (float *)calloc(numseeds, sizeof(float));
//     float * kseedsb = (float *)calloc(numseeds, sizeof(float));
//     float * kseedsx = (float *)calloc(numseeds, sizeof(float));
//     float * kseedsy = (float *)calloc(numseeds, sizeof(float));
//
//     GetLABXYSeeds_ForGivenStepSize(m_lvec, m_avec, m_bvec, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy,
//             STEP, m_height, m_width,xstrips,ystrips,  xerrperstrip, yerrperstrip,perturbseeds, edgemap);
//
//
//     PerformSuperpixelSLIC(m_lvec, m_avec, m_bvec, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP, edgemap,
//             compactness, numseeds, m_height, m_width, nb_iter, bres_lab, update_center, bres_contour, contour, offset);
//     numlabels = numseeds;
// //
//     int * nlabels = (int *)calloc(sz, sizeof(int));
//     EnforceLabelConnectivity(klabels, m_height, m_width, nlabels, &numlabels, (float)sz/((float)STEP*STEP), offset);
// //
//     for(int i = 0; i < sz; i++ )
//         klabels[i + offset] = nlabels[i];
//
//     //Free
//     free(m_lvec);
//     free(m_avec);
//     free(m_bvec);
//     free(nlabels);
//     free(kseedsl);
//     free(kseedsa);
//     free(kseedsb);
//     free(kseedsx);
//     free(kseedsy);
//     free(edgemap);
//
// }



void DoLABConversion01(
        float*      lvec,
        float*      avec,
        float*      bvec,
        float*      lvec01,
        float*      avec01,
        float*      bvec01,
        const int   m_height,
        const int   m_width)
{
    
    int sz = m_width*m_height;
    
    float max_l = -9999, max_a = -9999, max_b = -9999;
    float min_l = 9999, min_a = 9999, min_b = 9999;
    for( int j = 0; j < sz; j++ )
    {
        if (lvec[j] < min_l)
            min_l = lvec[j];
        if (avec[j] < min_a)
            min_a = avec[j];
        if (bvec[j] < min_b)
            min_b = bvec[j];
        if (lvec[j] > max_l)
            max_l = lvec[j];
        if (avec[j] > max_a)
            max_a = avec[j];
        if (bvec[j] > max_b)
            max_b = bvec[j];
    }
    for( int j = 0; j < sz; j++ )
    {
        lvec01[j] = max(min((lvec[j]-min_l)/(max_l-min_l),1),0);
        avec01[j] = max(min((avec[j]-min_a)/(max_a-min_a),1),0);
        bvec01[j] = max(min((bvec[j]-min_b)/(max_b-min_b),1),0);
    }
    
}

void DoSuperpixelSegmentation_ForGivenSuperpixelSize(
        float*                      ubuffr,
        float*                      ubuffg,
        float*                      ubuffb,
        const int					height,
        const int					width,
        int*						klabels,
        int                         numlabels,
        const int                   superpixelsize,
        const float                compactness,
        const int                   perturbseeds,
        int nb_iter,
        int bres_lab,
        int update_center,
        int bres_contour,
        float * contour,
        float * kernel,
        int pw,
        int h_bins)
{
    
    //------------------------------------------------
    const int STEP = sqrt((float)superpixelsize)+0.5;
    int m_width  = width;
    int m_height = height;
    int sz = m_width*m_height;
    
    
    float* m_lvec = (float *)calloc(sz, sizeof(float));
    float* m_avec = (float *)calloc(sz, sizeof(float));
    float* m_bvec = (float *)calloc(sz, sizeof(float));
    
    
    //--------------------------------------------------
    DoRGBtoLABConversion(ubuffr, ubuffg, ubuffb, m_lvec, m_avec, m_bvec, m_height, m_width);
    //--------------------------------------------------
    
    float* m_lvec01 = NULL;
    float* m_avec01 = NULL;
    float* m_bvec01 = NULL;
    
    if (bres_lab == 2) {
        m_lvec01 = (float *)calloc(sz, sizeof(float));
        m_avec01 = (float *)calloc(sz, sizeof(float));
        m_bvec01 = (float *)calloc(sz, sizeof(float));
        
        DoLABConversion01(m_lvec, m_avec, m_bvec, m_lvec01, m_avec01, m_bvec01, m_height, m_width);
    }
    
    // OR STAY IN RGB
//     for (int j=0; j<sz; j++) {
//         m_lvec[j]=ubuffr[j]*255;
//         m_avec[j]=ubuffg[j]*255;
//         m_bvec[j]=ubuffb[j]*255;
//     }
//
    int xstrips = 0.5+(float)m_width/(float)STEP;
    int ystrips = 0.5+(float)m_height/(float)STEP;
    int numseeds = xstrips*ystrips;
    
    int xerr = m_width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = m_width - STEP*xstrips;}
    int yerr = m_height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = m_height- STEP*ystrips;}
    
    float xerrperstrip = (float)xerr/(float)xstrips;
    float yerrperstrip = (float)yerr/(float)ystrips;
    
    float * kseedsl = (float *)calloc(numseeds, sizeof(float));
    float * kseedsa = (float *)calloc(numseeds, sizeof(float));
    float * kseedsb = (float *)calloc(numseeds, sizeof(float));
    float * kseedsx = (float *)calloc(numseeds, sizeof(float));
    float * kseedsy = (float *)calloc(numseeds, sizeof(float));
    
    GetLABXYSeeds_ForGivenStepSize(m_lvec, m_avec, m_bvec, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy,
            STEP, m_height, m_width,xstrips,ystrips,  xerrperstrip, yerrperstrip,perturbseeds);
    
//     if (compactness == 0) {
//         SLICO(m_lvec, m_avec, m_bvec, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP,
//                 compactness, numseeds, m_height, m_width, nb_iter, bres_lab, update_center, bres_contour, contour);
//     }
//     else {
    
    // precomputation patch
    int pos;
    int pos_i;
    float * Ip2 = (float *)calloc(sz*3, sizeof(float));
    float * Ip = (float *)calloc(sz*3, sizeof(float));
    for (int x=0; x<m_width; x++) {
        for (int y=0; y<m_height; y++) {
            float dist = 0;
            float count = 0;
            pos_i = y + x*m_height;
            for (int dx=-pw; dx<=pw; dx++) {
                for (int dy=-pw; dy<=pw; dy++) {
                    if ((x+dx<m_width)&&(x+dx>=0)&(y+dy>=0)&(y+dy<m_height)){
                        pos = y+dy + (x+dx)*m_height;
                        //float kk=kernel[dy+pw+(dx+pw)*(2*pw+1)];
                        float kk = 1;
                        Ip2[pos_i] += m_lvec[pos]*m_lvec[pos]*kk;
                        Ip2[pos_i+sz] += m_avec[pos]*m_avec[pos]*kk;
                        Ip2[pos_i+sz*2] += m_bvec[pos]*m_bvec[pos]*kk;
                        Ip[pos_i] += m_lvec[pos]*kk;
                        Ip[pos_i+sz] += m_avec[pos]*kk;
                        Ip[pos_i+sz*2] += m_bvec[pos]*kk;
                        count += kk;
                    }
                }
            }
            Ip2[pos_i] /= count;
            Ip2[pos_i+sz] /= count;
            Ip2[pos_i+sz*2] /= count;
            Ip[pos_i] /= count;
            Ip[pos_i+sz] /= count;
            Ip[pos_i+sz*2] /= count;
        }
    }
    
    
    SLIC(m_lvec, m_avec, m_bvec, kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP,
            compactness, numseeds, m_height, m_width, nb_iter, bres_lab, update_center, bres_contour, contour, pw,
            Ip2, Ip, h_bins, m_lvec01, m_avec01, m_bvec01);
//     }
    numlabels = numseeds;
//
    
    int * nlabels = (int *)calloc(sz, sizeof(int));
    EnforceLabelConnectivity(klabels, m_height, m_width, nlabels, &numlabels, (float)sz/((float)STEP*STEP));
// 
    for(int i = 0; i < sz; i++ )
        klabels[i] = nlabels[i];
    
    //Free
    free(Ip2);
    free(Ip);
    free(m_lvec);
    free(m_avec);
    free(m_bvec);
    free(m_lvec01);
    free(m_avec01);
    free(m_bvec01);
    free(nlabels);
    free(kseedsl);
    free(kseedsa);
    free(kseedsb);
    free(kseedsx);
    free(kseedsy);
    
}


void mexFunction(int nlhs, mxArray *plhs[],int nrhs,const mxArray *prhs[]){
    (void)nlhs;
    (void)nrhs;
    

    
    /* INPUTS */
    //RGB image
    float* img = (float*) mxGetPr(prhs[0]);
    
    int idx = 1;
    int SP_nbr = (int) mxGetScalar(prhs[idx++]);
    
    float compactness_i = (float) mxGetScalar(prhs[idx++]);
    int perturbseeds = 0;
    int nb_iter = 5;
    int bres_lab = 1;
    int update_center = 0;
    
    float compactness = (float) compactness_i;
    
    float * kernel = 0;
    
    float * contour; // = NULL;
    int bres_contour = 0;
    if (nrhs == 4) {
        bres_contour = 1;
        contour = (float*) mxGetPr(prhs[idx++]);
    }
    
    int pw = 0;
    int h_bins = 10;
    
    
    //Dimensions stuff
    const int* img_dims = mxGetDimensions(prhs[0]);
    int h_img = img_dims[0];
    int w_img = img_dims[1];
    int h_w_img = h_img*w_img;
    
    
    //Outputs
    int dims[2];
    dims[0] = h_img;
    dims[1] = w_img;
    plhs[0] = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    int * klabels = (int *)mxGetPr(plhs[0]);
    
    
    /*Tables*/
    //int *klabels = (int *)calloc(h_w_img, sizeof(int));
    float *ubuffr =(float *)calloc(h_w_img, sizeof(float));
    float *ubuffg =(float *)calloc(h_w_img, sizeof(float));
    float *ubuffb =(float *)calloc(h_w_img, sizeof(float));
    
    for(int i=0; i < h_w_img; i++) {
        ubuffr[i] = img[i];
        ubuffg[i] = img[i+h_w_img];
        ubuffb[i] = img[i+2*h_w_img];
        klabels[i] = -1;
    }
    
    
    const int superpixelsize = 0.5+(float)w_img*h_img/(float)SP_nbr;
    DoSuperpixelSegmentation_ForGivenSuperpixelSize(ubuffr, ubuffg, ubuffb, h_img, w_img ,
            klabels,SP_nbr,superpixelsize,compactness, perturbseeds, nb_iter, bres_lab, update_center, bres_contour, contour, kernel, pw, h_bins);
    
    
    free(ubuffr);
    free(ubuffg);
    free(ubuffb);
    
    
}







