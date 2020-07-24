/******************************************************************************
 * logging.c
 * -> Logging code for the QuickMAN SSE/SSE2-based Mandelbrot Set calculator
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
 * Initialize the logging structure.
 */
int log_init(logging* log)
{
	sprintf(log->file, LOG_FILE);
	return log_read(log, 1);
}

/**
 * Read a set of mandelbrot parms into a log entry (if not NULL).
 */
int log_read_entry(FILE* fp, log_entry* entry)
{
	double vals[5];
	int i = 0, j, pos;
	unsigned char strs[5][256], * str, c;

	// Read re, im, mag, iters, pal, and optional commands. To support legacy 
	// logfiles, the five main fields don't need any leading items (they can be 
	// just numbers).

	while (i < 5)
	{
		if (feof(fp))
			return 0;
		if (fgets((char*)&strs[i][0], sizeof(strs[0]), fp) == NULL)
			return 0;

		str = &strs[i][0];

		// Skip any leading whitespace
		pos = -1;
		do {
			c = str[++pos];
		} while (c == ' ' || c == '\t'); // null will terminate loop

		// ignore comments and new line symbol
		if (c == '/' || c == '#' || c == '\n')
			continue;

		// For added robustness, resync on "real", so corrupted files won't get 
		// us out of sync.
		if (!strncmp(&str[pos], "Real", 4))
			i = 0;

		// Might have an image parameter here (a number with or without leading 
		// items). Strip out any leading non-numeric/non-quote chars.
		for (j = pos; j < sizeof(strs[0]); j++)
		{
			c = str[j];
			if ((c >= '0' && c <= '9') || c == '-' || c == '.' || c == '\"' || !c)
				break;
		}
		if (c)
		{
			// Got something that looks like a number or a " if we get here. 
			// Any bad values will be set to 0.0 (ok). J is long lived (see 
			// below)
			vals[i] = atof(&str[j]);
			i++; // look for next value
		}
	}

	// All values good: update mandelbrot parms
	if (entry != NULL)
	{
		// Fill in the entry, including optional fields (if they're still at 
		// -1, nothing will happen later).
		entry->re = vals[0];
		entry->im = vals[1];
		entry->mag = vals[2];
		entry->max_iters = (unsigned int)vals[3];
		entry->palette = (unsigned int)vals[4];

		memset(entry->pal_name, 0, 256);

		// For user palette files (palette starts with " in logfile), use the 
		// position in the dropdown list. Assumes dropdown list is already 
		// populated. As a side effect this also allows the user to specify a 
		// builtin palette by either name (e.g. "Muted") or number (3)
		if (str[j] == '\"')
		{
			// Replace any trailing " with a null
			for (i = j + 1; i < sizeof(strs[0]); i++)
			{
				if (str[i] == '\"')
					str[i] = 0;
				if (!str[i])
					break;
			}
			sprintf(entry->pal_name, &str[j + 1]);
		}
	}
	return 1;
}

/**
 * Scan the logfile, dynamically allocate an array for the entries, and fill it 
 * in from the logfile. If init_pos is nonzero, initializes the position to the 
 * beginning.
 */
int log_read(logging* log, int init_pos)
{
	int i;
	FILE* fp;

	//_RPTN(_CRT_WARN, "log_read [name: %s, pos: %d]\n", log->file, init_pos);

	if (init_pos)
		log->pos = -1;

	// Kind of inefficient: scan once to get length, then scan again to fill in 
	// array.
	if ((fp = fopen(log->file, "rt")) == NULL)
		return 0;

	for (log->count = 0; log_read_entry(fp, NULL); log->count++);

	fclose(fp);

	if (!log->count)
		return 0; // normal cfg files will return here

	if (log->entries != NULL)   // Allocate the array and fill it in
		free(log->entries);

	if ((log->entries = (log_entry*)malloc(log->count * sizeof(log_entry))) == NULL)
		return 0;

	if ((fp = fopen(log->file, "rt")) == NULL)
		return 0;

	for (i = 0; i < log->count; i++)
		log_read_entry(fp, &log->entries[i]);

	fclose(fp);
	return 1;
}

/**
 * Open the logfile for appending and add the current image. Reset position if
 * reset_pos is 1.
 */
int log_update(man_calc_struct* m, logging* log, int reset_pos)
{
	char s[512], p[256];
	FILE* fp;

	if ((fp = fopen(log->file, "at")) == NULL) // open for append
		return 0;

	// For palette, use either number (for builtin palette), or "file" for user 
	// file.
	if (m->num_palettes > m->palette)
		sprintf_s(p, sizeof(p), "%d", m->palette);
	else
		sprintf_s(p, sizeof(p), "\"%s\"", m->pal_file);

	if (m->pal_xor) // add palette modification if it's in effect - v1.07
	{
		sprintf_s(s, sizeof(s), "\npal_xor 0x%06X", m->pal_xor);
		fputs(s, fp);
	}
	// Logfile read function ignores any leading items
	sprintf_s(s, sizeof(s),
		"\nReal    %-16.16lf"
		"\nImag    %-16.16lf"
		"\nMag     %-16lf"
		"\nIters   %d"
		"\nPalette %s\n",
		m->re, m->im, m->mag, m->max_iters, p);

	fputs(s, fp);
	fclose(fp);

	// Now reread the logfile (need to reallocate array) - a bit inefficient 
	// but who cares...
	// Keep current position
	return log_read(log, reset_pos);
}

/**
 * Get the next or prev entry from the log entry array. Returns the entry.
 */
log_entry* get_log_entry(man_calc_struct* m, logging* log, int next_prevn)
{
	log_entry* entry;

	if (log->entries == NULL)
		return NULL;

	if (next_prevn)
	{
		if (++log->pos > log->count - 1)
			// log_pos = log_count - 1;	// stop at end
			log->pos = 0;				// wrap to beginning
	}
	else if (--log->pos < 0)
			// log_pos = 0;				// stop at beginning
			log->pos = log->count - 1;	// wrap to end

	entry = &log->entries[log->pos];

	m->re = entry->re;
	m->im = entry->im;
	m->mag = entry->mag;
	m->max_iters = entry->max_iters;
	if (!(m->status & STAT_PALETTE_LOCKED))
		m->palette = entry->palette;

	return entry;
}

/**
 * Set palette number.
 *
 * unused
 */
void set_palette_number(logging* log, char* pal_name, int palette)
{
	int i;
	size_t sz;

	if (log->entries == NULL)
		return;

	for (i = 0; i < log->count; i++)
	{
		if (log->entries[i].pal_name[0] == 0)
			continue;

		sz = strlen(pal_name);
		if (strncmp(pal_name, log->entries[i].pal_name, sz) == 0)
		{
			log->entries[i].palette = palette;
			break;
		}
	}
}