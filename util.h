/*
 * util.h: Utilitary functions.
 *
 * Copyright (C) 2012 Mikhail Kurnosov
 */

#ifndef UTIL_H
#define UTIL_H

void *xmalloc(size_t size);
void *xrealloc(void *ptr, size_t size);
void exit_error(const char *format, ...);

#endif /* UTIL_H */
