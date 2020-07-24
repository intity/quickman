/******************************************************************************
 * cfg.c
 * -> CFG code for the QuickMAN SSE/SSE2-based Mandelbrot Set calculator
 * Copyright (C) 2006-2008 Paul Gentieu (paul.gentieu@yahoo.com)
 *
 * This file is part of QuickMAN.
 *
 * QuickMAN is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
 *
 * Project Page: http://quickman.sourceforge.net
 *
 * Author: Paul Gentieu (main code, ASM iteration cores, palettes)
 *
 ******************************************************************************
 */

#define STRICT

#include "quickman.h"
#include <stdlib.h>
#include <string.h>

/**
 * Initializes the configuration file.
 * Call this early, before the main window is initialized could change window
 * size.
 */
int cfg_init(man_calc_struct* m, settings* cfg, char* file)
{
	FILE* fp;
	int i, pos, val;
	size_t sz;
	char str[256], c;
	setting* s = (setting*)cfg; // treat structs as arrays

	if ((fp = fopen(file, "rt")) == NULL)
		return 0;

	while (1)
	{
		if (feof(fp))
			break;
		if (fgets(str, 256, fp) == NULL)
			break;

		// Skip any leading whitespace
		pos = -1;
		do {
			c = str[++pos];
		} while (c == ' ' || c == '\t'); // null will terminate loop

		// ignore comments and new line symbol
		if (c == '/' || c == '#' || c == '\n')
			continue;

		for (i = 0; i < CFG_SIZE; i++)
		{
			sz = strlen(s[i].name);
			if (strncmp(s[i].name, str, sz) == 0)
			{
				// Can use the function for reading a palette RGB value here 
				// (it will read normal integers too). As a side effect any 24 
				// bit value can be specified as three individual bytes if 
				// desired...
				// Some settings are palette values.

				get_palette_rgb_val(pos + (int)sz, str, 256, &val);

				// set value if it's within legal range
				if (val >= s[i].min && val <= s[i].max)
					s[i].c_val = val;
				break;
			}
		}
	}
	fclose(fp);

	// maybe eliminate these separate variables later
	m->xsize = m->prev_xsize = cfg->xsize.c_val;
	m->ysize = m->prev_ysize = cfg->ysize.c_val;

	// just so 1st dialog box doesn't get div by 0
	m->min_dimension = m->xsize > m->ysize ? m->ysize : m->xsize;
	m->max_iters_color = cfg->max_iters_color.c_val;

	return 1;
}

/**
 * Autoreset settings fields to defaults, if so configured. Call before every
 * image recalculation. Only should be called with global cfg_settings (change
 * to not take parm?)
 *
 * unused
 */
void autoreset_settings(settings* cfg)
{
	int i;
	setting* d;
	d = (setting*)cfg; // treat struct as array
	for (i = 0; i < CFG_SIZE; i++)
		if (SETTING_AUTORESET(&d[i]))
			d[i].c_val = d[i].d_val;
}

/**
 * Set all fields to -1 value, which invalidates all settings.
 *
 * unused
 */
void invalidate_settings(settings* cfg)
{
	int i;
	setting* d;

	d = (setting*)cfg; // treat struct as array
	for (i = 0; i < CFG_SIZE; i++)
		d[i].c_val = -1;
}

/**
 * Copy any settings fields that have changed (i.e., are >= 0) from src to dest.
 * If copy_to_default is 1, also copies changed settings to the default_val
 * fields (use with quickman.cfg). Then autoreset will reset to the
 * quickman.cfg default values.
 *
 * unused
 */
void copy_changed_settings(settings* dest, settings* src, int copy_to_default)
{
	int i;
	setting* s;
	setting* d;

	s = (setting*)src;	// treat structs as arrays
	d = (setting*)dest;
	for (i = 0; i < CFG_SIZE; i++)
		if (s[i].c_val >= 0)
		{
			d[i].c_val = s[i].c_val;
			if (copy_to_default)
				d[i].d_val = s[i].c_val;
		}
}