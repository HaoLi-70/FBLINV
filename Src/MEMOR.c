
#include "MEMOR.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:
        8 Sept. 2021.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

static int inster_scalar(STR_SCALAR *ptr, double val);

static int inster_1d(STR_1D *ptr, double val, int indx1);

static int inster_2d(STR_2D *ptr, double val, int indx1, int indx2);

static int inster_3d(STR_3D *ptr, double val, int indx1, int indx2, int indx3);

static int inster_4d(STR_4D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4);

static int inster_5d(STR_5D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5);

static int inster_7d(STR_7D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6, int indx7);

static int inster_8d(STR_8D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6, int indx7, int indx8);

static double *ele1d(STR_1D *ptr, int indx1);

static double *ele2d(STR_2D *ptr, int indx1, int indx2);

static double *ele3d(STR_3D *ptr, int indx1, int indx2, int indx3);

static double *ele4d(STR_4D *ptr, int indx1, int indx2, int indx3, int indx4);

static double *ele5d(STR_5D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5);

static double *ele7d(STR_7D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6, int indx7);

static double *ele8d(STR_8D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6, int indx7, int indx8);

/*----------------------------------------------------------------------------*/

static int inster_scalar(STR_SCALAR *ptr, double val){
    
    if(ptr->dat == NULL){
      ptr->dat = (double *)malloc(sizeof(double));

      ptr->dat[0] = val;
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static int inster_1d(STR_1D *ptr, double val, int indx1){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_SCALAR *)malloc(sizeof(STR_SCALAR) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_scalar(ptr->dat+indx1-ptr->bounds[0], val);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_scalar(ptr->dat+indx1-ptr->bounds[0], val);

      }else{

        STR_SCALAR *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_SCALAR *)malloc(sizeof(STR_SCALAR) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_scalar(ptr->dat+indx1-ptr->bounds[0], val);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static int inster_2d(STR_2D *ptr, double val, int indx1, int indx2){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_1D *)malloc(sizeof(STR_1D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_1d(ptr->dat+indx1-ptr->bounds[0], val, indx2);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_1d(ptr->dat+indx1-ptr->bounds[0], val, indx2);

      }else{

        STR_1D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_1D *)malloc(sizeof(STR_1D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_1d(ptr->dat+indx1-ptr->bounds[0], val, indx2);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static int inster_3d(STR_3D *ptr, double val, int indx1, int indx2, int indx3){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_2D *)malloc(sizeof(STR_2D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_2d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_2d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3);

      }else{

        STR_2D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_2D *)malloc(sizeof(STR_2D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_2d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static int inster_4d(STR_4D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_3D *)malloc(sizeof(STR_3D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_3d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_3d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4);

      }else{

        STR_3D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_3D *)malloc(sizeof(STR_3D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_3d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static int inster_5d(STR_5D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5){


    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_4D *)malloc(sizeof(STR_4D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_4d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
          indx5);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_4d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5);

      }else{

        STR_4D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_4D *)malloc(sizeof(STR_4D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_4d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, indx5);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

extern int inster_6d(STR_6D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_5D *)malloc(sizeof(STR_5D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_5d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
          indx5, indx6);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_5d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6);

      }else{

        STR_5D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_5D *)malloc(sizeof(STR_5D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_5d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static int inster_7d(STR_7D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6, int indx7){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_6D *)malloc(sizeof(STR_6D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_6d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
          indx5, indx6, indx7);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_6d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6, indx7);

      }else{

        STR_6D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_6D *)malloc(sizeof(STR_6D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_6d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6, indx7);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static int inster_8d(STR_8D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6, int indx7, int indx8){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_7D *)malloc(sizeof(STR_7D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_7d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
          indx5, indx6, indx7, indx8);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_7d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6, indx7, indx8);

      }else{

        STR_7D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_7D *)malloc(sizeof(STR_7D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_7d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6, indx7, indx8);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

extern int inster_9d(STR_9D *ptr, double val, int indx1, int indx2, int indx3, \
    int indx4, int indx5, int indx6, int indx7, int indx8, int indx9){

    int i;

    if(ptr->dat == NULL){
      if(indx1<0){
        ptr->bounds[0] = indx1;
        ptr->bounds[1] = 0;

      }else{
        ptr->bounds[0] = 0;
        ptr->bounds[1] = indx1;
      }
      ptr->dat = (STR_8D *)malloc(sizeof(STR_8D) \
          *(ptr->bounds[1]-ptr->bounds[0]+1));

      for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
        ptr->dat[i].dat = NULL;
      }

      inster_8d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
          indx5, indx6, indx7, indx8, indx9);

    }else{
      if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){

        inster_8d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6, indx7, indx8, indx9);

      }else{

        STR_8D *ptr_tmp;
        int bounds[2];
        int tmp;

        if(indx1<ptr->bounds[0]){
          bounds[0] = indx1;
          bounds[1] = ptr->bounds[1];
          tmp = bounds[0]-indx1;
        }else{
          bounds[0] = ptr->bounds[0];
          bounds[1] = indx1;
          tmp = 0;
        }

        ptr_tmp = (STR_8D *)malloc(sizeof(STR_8D) \
            *(bounds[1]-bounds[0]+1));

        for(i=0;i<=bounds[1]-bounds[0]+1;i++){
          ptr_tmp[i].dat = NULL;
        }

        for(i=0;i<=ptr->bounds[1]-ptr->bounds[0]+1;i++){
          ptr_tmp[i+tmp].dat = ptr->dat[i].dat;
        }

        free(ptr->dat);
        ptr->dat = ptr_tmp;

        ptr->bounds[0] = bounds[0];
        ptr->bounds[1] = bounds[1];

        inster_8d(ptr->dat+indx1-ptr->bounds[0], val, indx2, indx3, indx4, \
            indx5, indx6, indx7, indx8, indx9);

      }
    }

    return 0;

}

/*----------------------------------------------------------------------------*/

static double *ele1d(STR_1D *ptr, int indx1){

    if(ptr->dat == NULL) return NULL;

//    fprintf(stderr," 1d %d %d %d \n",ptr->bounds[0],ptr->bounds[1], indx1);

    if(indx1 > ptr->bounds[0] && indx1 < ptr->bounds[1]){
      return NULL;
    }else{
      return ptr->dat[indx1].dat;
    }
}

/*----------------------------------------------------------------------------*/

static double *ele2d(STR_2D *ptr, int indx1, int indx2){

    if(ptr->dat == NULL) return NULL;

//    fprintf(stderr," 2d %d %d %d \n",ptr->bounds[0],ptr->bounds[1], indx1);

    if(indx1 > ptr->bounds[0] && indx1 < ptr->bounds[1]){
      return NULL;
    }else{
      return ele1d(ptr->dat+indx1-ptr->bounds[0], indx2);
    }
}

/*----------------------------------------------------------------------------*/

static double *ele3d(STR_3D *ptr, int indx1, int indx2, int indx3){

    if(ptr->dat == NULL) return NULL;

//    fprintf(stderr," 3d %d %d %d \n",ptr->bounds[0],ptr->bounds[1], indx1);

    if(indx1 > ptr->bounds[0] && indx1 < ptr->bounds[1]){
      return NULL;
    }else{
      return ele2d(ptr->dat+indx1-ptr->bounds[0], indx2, indx3);
    }
}

/*----------------------------------------------------------------------------*/

static double *ele4d(STR_4D *ptr, int indx1, int indx2, int indx3, int indx4){

    if(ptr->dat == NULL) return NULL;

//    fprintf(stderr," 4d %d %d %d \n",ptr->bounds[0],ptr->bounds[1], indx1);

    if(indx1 > ptr->bounds[0] && indx1 < ptr->bounds[1]){
      return NULL;
    }else{
      return ele3d(ptr->dat+indx1-ptr->bounds[0], indx2, indx3, indx4);
    }
}

/*----------------------------------------------------------------------------*/

static double *ele5d(STR_5D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5){

    if(ptr->dat == NULL) return NULL;

//    fprintf(stderr," 5d %d %d %d \n",ptr->bounds[0],ptr->bounds[1], indx1);

    if(indx1 > ptr->bounds[0] && indx1 < ptr->bounds[1]){
      return NULL;
    }else{
      return ele4d(ptr->dat+indx1-ptr->bounds[0], indx2, indx3, indx4, \
          indx5);
    }
}

/*----------------------------------------------------------------------------*/

extern double *ele6d(STR_6D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6){

    if(ptr->dat == NULL) return NULL;

//    fprintf(stderr," 6d %d %d %d \n",ptr->bounds[0],ptr->bounds[1], indx1);

    if(indx1 > ptr->bounds[0] && indx1 < ptr->bounds[1]){
      return NULL;
    }else{
      return ele5d(ptr->dat+indx1-ptr->bounds[0], indx2, indx3, indx4, \
          indx5, indx6);
    }
}

/*----------------------------------------------------------------------------*/

static double *ele7d(STR_7D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6, int indx7){

    if(ptr->dat == NULL) return NULL;

    if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){
      return NULL;
    }else{
      return ele6d(ptr->dat+indx1-ptr->bounds[0], indx2, indx3, indx4, \
          indx5, indx6, indx7);
    }
}

/*----------------------------------------------------------------------------*/

static double *ele8d(STR_8D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6, int indx7, int indx8){

    if(ptr->dat == NULL) return NULL;

    if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){
      return NULL;
    }else{
      return ele7d(ptr->dat+indx1-ptr->bounds[0], indx2, indx3, indx4, \
          indx5, indx6, indx7, indx8);
    }
}

/*----------------------------------------------------------------------------*/

extern double *ele9d(STR_9D *ptr, int indx1, int indx2, int indx3, int indx4, 
    int indx5, int indx6, int indx7, int indx8, int indx9){

    if(ptr->dat == NULL) return NULL;

    if(indx1 >= ptr->bounds[0] && indx1 <= ptr->bounds[1]){
      return NULL;
    }else{
      return ele8d(ptr->dat+indx1-ptr->bounds[0], indx2, indx3, indx4, \
          indx5, indx6, indx7, indx8, indx9);
    }
}

/*----------------------------------------------------------------------------*/
