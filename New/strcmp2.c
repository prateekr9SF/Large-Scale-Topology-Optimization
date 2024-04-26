/** 
 * @file strcmp2.c
 * @brief Implementation of a custom string comparison function for CalculiX.
 *
 * CalculiX - A 3-dimensional finite element program.
 * This file contains a custom implementation of a string comparison function,
 * which compares up to a specified number of characters of two strings, or less
 * if one of the strings ends before reaching that number. This is especially 
 * useful in finite element methods where precision and control over string
 * comparison lengths are required.
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
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "main.h"

/** 
 * @brief Compares two strings up to a specified number of characters.
 *
 * This function compares the string pointed to by s1 to the string pointed to by s2,
 * but only up to 'length' characters or until the end of either string is reached.
 * This custom behavior ensures that comparisons are made only within the bounds of
 * actual string content, avoiding overflows and mismatches caused by excess comparison.
 *
 * @param s1 Pointer to the first null-terminated string to be compared.
 * @param s2 Pointer to the second null-terminated string to be compared.
 * @param length The maximum number of characters to compare.
 * @return An integer less than, equal to, or greater than zero if `s1` is found,
 * respectively, to be less than, to match, or be greater than `s2`.
 */
ITG strcmp2(const char *s1, const char *s2, ITG length)
{
    ITG a, b, i = 0;

    do {
        a = *s1++;
        b = *s2++;

        // Terminate comparison if end of either string is reached
        if (b == '\0' || a == '\0') {
            a = '\0';
            b = '\0';
            break;
        }
        i++;
    } while (a == b && i < length);

    return (a - b);
}
