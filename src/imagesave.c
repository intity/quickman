/******************************************************************************
 * imagesave.c
 * -> Code for saving QuickMAN images.
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
 *
 * 11/09/08 PG: New file for v1.10.
 */

#define STRICT

#include <png.h>
#include <setjmp.h>

#include "quickman.h"

/* ------------------------------ PNG functions ---------------------------- */

static png_structp png;
static png_infop pnginfo;
static FILE* png_fp;

/**
 * Start the PNG save for an image of dimensions WIDTH x HEIGHT.
 */
int png_save_start(char* file, int width, int height)
{
	png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

	if (png == NULL)
	{
		png_error(NULL, "Failed to create png object.");
		return 0;
	}

	if ((pnginfo = png_create_info_struct(png)) == NULL)
	{
		png_error(png, "Failed to create png_info object.");
		png_destroy_write_struct(&png, NULL);
		return 0;
	}

	if (fopen_s(&png_fp, file, "wb"))
	{
		png_error(png, "Failed to open png file.");
		png_destroy_write_struct(&png, &pnginfo);
		return 0;
	}

	if (setjmp(png_jmpbuf(png)))
	{
		png_error(png, "Any pnglib errors jump here.");
		png_destroy_write_struct(&png, &pnginfo);
		fclose(png_fp);
		return 0;
	}

	png_init_io(png, png_fp);
	
	// Set B,G,R color order: necessary for correct palette mapping
	png_set_bgr(png);

	// Set Zlib compression levels; 0 = none; 9 = best. 3-6 are supposed to be 
	// almost as good as 9 for images. Test: 6- 43s 17519kb; 9- 45s 18573kb; 9 
	// sometimes gives worse compression
	png_set_compression_level(png, 6);

	png_set_IHDR(png, pnginfo, width, height, 8, 
		PNG_COLOR_TYPE_RGB, 
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT, 
		PNG_FILTER_TYPE_DEFAULT);

	png_write_info(png, pnginfo);
	return 1;
}

/**
 * Write one row of the image. Row must have 3 * width bytes (RGB), where width
 * was the parameter passed to png_save_start.
 */
int png_save_write_row(unsigned char* row)
{
	if (setjmp(png_jmpbuf(png)))
	{
		png_error(png, "Any pnglib errors jump here.");
		png_destroy_write_struct(&png, &pnginfo);
		fclose(png_fp);
		return 0;
	}
	png_write_row(png, row);
	return 1;
}

/**
 * Finish the PNG save. Call after all rows are written.
 */
int png_save_end(void)
{
	if (setjmp(png_jmpbuf(png)))
	{
		png_error(png, "Any pnglib errors jump here.");
		png_destroy_write_struct(&png, &pnginfo);
		fclose(png_fp);
		return 0;
	}

	png_write_end(png, NULL);
	png_destroy_write_struct(&png, &pnginfo);
	fclose(png_fp);
	return 1;
}

// Later add functions for saving raw RGB data, magnitudes, and iteration 
// counts
