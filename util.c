/*
 * util.c: Utilitary functions.
 *
 * Copyright (C) 2012 Mikhail Kurnosov
 */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>

#include "util.h"

enum {
    ERROR_BUF_MAX = 512
};

/* xmalloc: */
void *xmalloc(size_t size)
{
    void *p;

    if (size == 0)
        return NULL;

    if ( (p = malloc(size)) == NULL) {
        exit_error("xmalloc: No enough memory");
    }
    return p;
}

/* xrealloc: */
void *xrealloc(void *ptr, size_t size)
{
    void *p;

    if ( (p = realloc(ptr, size)) == NULL) {
        exit_error("xrealloc: No enough memory");
    }
    return p;
}

/* exit_error: Prints an error message and terminates. */
void exit_error(const char *format, ...)
{
    va_list ap;
    static char buf[ERROR_BUF_MAX];

    va_start(ap, format);
    vsprintf(buf, format, ap);
    va_end(ap);
    if (strlen(buf) > 0)
        fprintf(stderr, "error: %s\n", buf);
    exit(EXIT_FAILURE);
}

