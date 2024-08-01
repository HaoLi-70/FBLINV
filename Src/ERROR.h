
#ifndef ERROR_h
#define ERROR_h

/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

/*----------------------------------------------------------------------------*/

enum error_level {enum_warning, enum_error};

#define MAX_MESSAGE_LENGTH 2000

extern char Err_Mes[MAX_MESSAGE_LENGTH];

/*----------------------------------------------------------------------------*/

extern int Error(enum error_level error_lv, const char *routine_name,
    const char *message_str);

extern bool FILE_EXIST(char *Filename);

/*----------------------------------------------------------------------------*/

#endif /* ERROR_h */
