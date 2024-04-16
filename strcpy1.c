/** 
 * @file main.c
 * @brief Implementation of custom string operations for CalculiX.
 *
 * CalculiX - A 3-dimensional finite element program.
 * This file contains implementations of specific string operations 
 * that are tailored to the needs of the CalculiX program, including
 * a customized string copy function that pads with spaces if the source
 * string is shorter than the desired length.
 *
 * Copyright (C) 1998-2018 Guido Dhondt.
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation (version 2).
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "main.h"

/** 
 * @brief Copies a string from source to destination with space padding.
 *
 * This function copies a string from 's2' to 's1' up to 'length' characters.
 * If 's2' is shorter than 'length', the remainder of 's1' is filled with spaces.
 * This is particularly useful in scenarios where strings need to be a fixed width.
 *
 * @param s1 Pointer to the destination string buffer.
 * @param s2 Pointer to the source string to be copied.
 * @param length The number of characters to copy or pad.
 * @return Always returns 0.
 */
ITG strcpy1(char *s1, const char *s2, ITG length)
{
    ITG b, i, blank = 0;

    for (i = 0; i < length; i++) {
        if (blank == 0) {
            b = *s2;
            if (b == '\0') blank = 1;
        }
        if (blank == 0) { *s1 = *s2; s2++; }
        else *s1 = ' ';
        s1++;
    }
    return 0;
}




