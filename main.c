/*
 * main.c: Main module.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov
 */

#include "md.h"

int main(int argc, char **argv)
{
    if (argc == 3) {
        md_initialize(argv[1], argv[2]);
    } else if (argc == 2) {
        md_initialize(argv[1], "md.input");
    } else {
        md_initialize("mdf.cls", "md.input");
    }
    md_run();
    md_finalize();

    return 0;
}
