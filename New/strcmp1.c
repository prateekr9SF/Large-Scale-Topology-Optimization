/** 
 * @file strcmp1.c
 * @brief Implementation of a custom string comparison function for CalculiX.
 *
 * CalculiX - A 3-dimensional finite element program.
 * This file contains a custom string comparison function which
 * handles comparisons particularly suited for scenarios where one string
 * is variable-length and the other is a fixed-length, as typically found
 * in data parsing scenarios in finite element analyses.
 *
 * Copyright (C) 1998-2018 Guido Dhondt.
 * This program is licensed under the GNU General Public License version 2.
 * It is distributed without any warranty; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
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
 * @brief Compares two strings with special handling for variable and fixed-length strings.
 *
 * Compares the string pointed to by s1 to the string pointed to by s2.
 * Unlike standard strcmp, this function terminates the comparison if the end
 * of either string (variable or fixed) is reached, ensuring that s1 does not
 * exceed its intended length.
 *
 * @param s1 Pointer to the first null-terminated string to be compared.
 * @param s2 Pointer to the second null-terminated string to be compared.
 * @return An integer less than, equal to, or greater than zero if `s1` is found,
 * respectively, to be less than, to match, or be greater than `s2`.
 */
ITG strcmp1(const char *s1, const char *s2)
{
    ITG a, b;

    do {
        a = *s1++;
        b = *s2++;

        // Handle cases where one or both strings reach their end.
        if (b == '\0') {
            a = '\0';
            b = '\0';
            break;
        }
        if (a == '\0') {
            a = '\0';
            b = '\0';
            break;
        }
    } while (a == b);

    return (a - b);
}
