
#include "ERROR.h"

/*----------------------------------------------------------------------------*/

    /*######################################################################
     
     revision log:
        8 Sept. 2021.
     
    ######################################################################*/

/*----------------------------------------------------------------------------*/

char ERR_MESS[MAX_MESSAGE_LENGTH] = "";

extern int Error(enum error_level error_lv, const char *routine_name,\
    const char *message_str){
    
    /*######################################################################
      Purpose:
        error handler.
      Record of revisions:
        16 Jan. 2023.
      Input parameters:
        error_lv, level of the error (enum_warning, enum_error)
        routine_name, name of the routine.
        message_str, the message.
      Return:
        .
    ######################################################################*/

    switch (error_lv) { 
      case enum_error:
        fprintf(stderr, "\n\n-ERROR in routine %s\n %s\n", \
            routine_name, (message_str) ? message_str : " (Undocumented)\n\n");
        MPI_Finalize();
        exit(0);
        return 0;
      case enum_warning:
        fprintf(stderr, "\n\n-WARNING in routine %s\n %s\n", routine_name, \
            (message_str) ? message_str : " (Undocumented)\n\n");
          return 0;
      default:
        return 1;
    }
    return 0;
}

/*----------------------------------------------------------------------------*/

extern bool FILE_EXIST(char *Filename){
      
    /*######################################################################
      Purpose:
        check if a file exists.
      Record of revisions:
        23 Aug. 2023.
      Input parameters:
        Filename, name of the file
      Return:
        true or false
    ######################################################################*/

    FILE *Fa = fopen(Filename, "r");

    if(Fa){

      fclose(Fa);

      return true;
    }

    return false;
}

/*----------------------------------------------------------------------------*/

