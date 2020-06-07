/******************************************************************************
 * quickman.c
 * -> Main file for the QuickMAN SSE/SSE2-based Mandelbrot Set calculator
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
 * 11/05/08 - 11/18/08 PG: v1.10 (started from v1.07)
 *
 * Removed previous history comments. See older versions for those.
 *
 * Added Save Image functionality and dialog box controls. Currently can save 
 * arbitrary-sized PNG images. Eventually will be able to save iteration and 
 * magnitude data as well (the selection checkboxes are currently grayed out).
 *
 * Uses a stripped-down version of libpng to implement the PNG save 
 * functionality. A compiled library and build instructions is included in the 
 * QuickMAN distribution.
 *
 * Modified the calculation and rendering engines to run from a control 
 * structure (man_calc_struct), so multiple instances can run in separate 
 * threads. This allows the save function to run in the background while normal 
 * browsing continues. See more comments at do_save().
 *
 * Removed some obsolete comments.
 * Removed Win98 support.
 * Split off most structure definitions and defines into the quickman.h header 
 * file.
 * Reduced default pan rate to better suit current CPUs.
 * Added quickman.cfg options bits to allow normalized rendering and the exact 
 * algorithm to be chosen as the defaults.
 *
 * Bugs and annoyances fixed:
 *
 * The max_iters color was being inverted with the rest of the palette. Wasn't 
 * the intention fixed Wasn't displaying "Logged" in the status line after the 
 * log image button was pressed.
 *
 * TODO: Fix palette locking/inverted indicator can be confusing
 */

#define STRICT

#include "quickman.h"

#include <stdlib.h>
#include <mmsystem.h>	// for timer functions
#include <process.h>	// for threading functions
#include <math.h>

#ifdef _WIN32
#include <commctrl.h>
	#ifdef _DEBUG
	#include <crtdbg.h>
	#endif
#endif

/**
 * Combo box initialization strings/defines. String order should correspond
 * to the above define order for PRECISION_* and ALG_*.
 */
static char* precision_strs[] = { "Auto", "Single", "Double", "Extended" };
static char* alg_strs[] = { 
	"Fast, AMD",   "Exact, AMD",
	"Fast, Intel", "Exact, Intel",
	"Fast, C",     "Exact, C" 
};
//	"Fast error" };

/**
 * Striped, Flaming+, and Plantlike are marked for replacement- rarely used.
 * Removed leading numbers to give more space for palette names.
 */
static char* palette_strs[] = {
	"Monochrome", 
	"Striped", 
	"Loud", 
	"Muted", 
	"Purple",
	"Earthy", 
	"Smoky", 
	"Acid", 
	"Flaming", 
	"Metallic",
	"Angry", 
	"Dreamy", 
	"Flaming+", 
	"Plantlike" 
};

static char* rendering_strs[] = { "Standard", "Normalized" };

// Not all of these will be used
static char* num_threads_strs[] = { 
	"1", "2", "4", "8", "16", "32", "64", "128", "256" 
};

static const char help_text[] = {
	"For complete documentation, please go to the QuickMAN project webpage\n"
	"and click on the Documentation tab.\n\n"

	"http://quickman.sourceforge.net\n\n"

	"Operation Summary\n\n"
	"\tMouse buttons: zoom in/out; zoom rectangle in magnifier mode\n"
	"\tMouse wheel: increase/decrease Max Iters\n\n"

	"Keybindings\n\n"
	"F1\t: show this message\n"
	"Arrow keys or A, S, D, W: move around the image (pan)\n"
	"Space\t: (with mouse): drag the image\n"
	"Shift\t: (with/without arrow keys): start/stop automatic panning\n"
	"Ctrl\t: (during panning): increase panning speed\n"
	"F\t: or [Fullscreen] button: switch between windowed and fullscreen\n"
	"Esc\t: exit fullscreen mode\n"
	"C\t: show/hide the control window\n"
	"N\t: or Next button: go to the next logfile image\n"
	"P\t: or Previous button: go to the previous logfile image\n"
	"H\t: or [Home] button: go to the home image\n"
	"L\t: lock the current palette (ignore logfile palettes)\n"
	"I\t: invert the current palette\n"
	"Z\t: switch between realtime zooming and magnifier modes\n"
};

// Global mandelbrot parameters
static int num_threads_ind = 0;		// number of threads string index
int num_threads = 1;				// number of calculation threads. 
									// Limited to values that make sense: 
									// 1, 2, 4, 8, 16...
static int prev_xsize;				// previous sizes, for restoring window
static int prev_ysize;
static double mouse_re;				// re/im coordinates of the mouse position
static double mouse_im;
static double zoom_start_mag;		// starting magnification, for zoom button
static int precision_loss = 0;			// 1 if precision loss detected on most 
										// recent calculation
static unsigned int num_builtin_palettes;	// number of builtin palettes
static unsigned int num_palettes;			// total number of palettes
static char palette_file[256];		// filename of current user palette file
static char logfile[256] = LOG_FILE;	// default logfile
static char savefile[256] = "image_1";	// default save filename

/**
 * Presets for the logfile combo box. Only need to add files here that we want 
 * to be at the top of the box. Others will be added from the current directory.
 */
char* file_strs[] = { logfile, "auto_panzoom.log" };

// Timing/benchmarking values
static double iter_time = 0.0;          // mandelbrot iteration time only
static double zoom_time;                // time to do 1 zoom, includes above + 
										// extra overhead
static double calc_interval_time = 0.0; // mandelbrot time during one interval
static double calc_total_time = 0.0;    // mandelbrot total time
static double interval_time = 0.0;      // pan/zoom time during 1 interval
static double total_time = 0.0;         // total pan/zoom time
static unsigned int total_frames = 0;       // total frames
static unsigned int interval_frames = 0;    // frames per interval
static double file_tot_time = 0.0;      // total calculation time since file 
										// opened
static int all_recalculated = 0;        // flag indicating a recalculation of 
										// the whole image happened
static TIME_UNIT zoom_start_time;       // for zoom button benchmarking (helps 
										// measure overhead)

// Idicates whether processor supports SSE/SSE2 and CMOV instructions
static int sse_support = 0; // 1 for SSE, 2 for SSE and SSE2

// Mouse position: index 0 is position on initial button press; index 1 is 
// current position
static int mouse_x[2], mouse_y[2];

static int nav_mode = MODE_RTZOOM;
static int do_rtzoom = 0;		// if nonzero, do a realtime zoom in/out
static int prev_do_rtzoom = 0;	// previous state of do_rtzoom

/**
 * Zoom steps as a function of slider value. Nonlinear at beginning and end.
 * These seem to work pretty well.
 */
static double rtzoom_mag_steps[] = {
	1.000625, 	// 0, super slow
	1.001250, 
	1.002500, 
	1.005000, 
	1.010000, 
	1.015000, 
	1.020000, 
	1.025000, 
	1.030000, 
	1.040000, 
	1.050000, 	// 10, default. 1.05 was previous value, but no existing 
				// benchmarks, so can change it
	1.060000, 
	1.070000, 
	1.080000, 
	1.090000, 
	1.100000, 
	1.110000, 
	1.120000, 
	1.140000, 
	1.170000,
	1.200000	// 20, fast
};

// Arbitrary units; slider range is determined by array size
#define MAX_ZOOM_RATE	(NUM_ELEM(rtzoom_mag_steps) - 1)
#define DEF_ZOOM_RATE	(MAX_ZOOM_RATE >> 1) // default zoom rate

/**
 * Pan step scales as a function of slider value.
 */
static double pan_step_scales[] = {
   0.001250,	// 0, very slow
   0.002500, 
   0.005000, 
   0.010000, 
   0.020000, 
   0.040000, 
   0.080000, 
   0.200000, 
   0.400000, 
   0.600000,
   0.800000,	// 10, default. With new magic constants, this should be 
				// roughly compatible with prev. benchmarks
   1.000000, 
   1.200000, 
   1.400000, 
   1.600000, 
   1.800000, 
   2.000000, 
   2.200000, 
   2.400000, 
   2.600000,
   2.800000		// 20, fast
};

#define MAX_PAN_RATE	(NUM_ELEM(pan_step_scales) - 1)
#define DEF_PAN_RATE	(MAX_PAN_RATE >> 1) // default pan rate

// GUI stuff
static HCURSOR mag_cursor, rtzoom_cursor, hopen_cursor, hclosed_cursor;
static HCURSOR arrow_cursor, wait_cursor;
static HCURSOR mag_zoom_cursor;		// either mag or rtzoom, depending on mode
static RECT main_rect;				// rectangle for the main window
static HWND hwnd_main;				// handle for the main window
static HWND hwnd_dialog = 0;		// handle for the dialog box
static HWND hwnd_info, hwnd_status, hwnd_status2;	// handles for the 
													// info/status text areas
static HWND hwnd_iters, hwnd_thumbnail_frame;		// and other controls
static HINSTANCE hinstance = 0;		// the application
static HDC hscreen_dc = NULL;		// screen device context
static int x_border;		// variables for calculating window sizes and such
static int y_border;
static int y_thinborder;
static int x_dialog_border;
static int y_dialog_border;
static int lpix_per_inch;

static int status = 0;		// general status bitfield, sstat

/**
 * Quadrant-based panning algorithm:
 *
 * Each quadrant is the size of the screen (xsize x ysize). Initially, the 
 * screen window is placed in the upper left (UL) quadrant, at position (0, 0).
 *
 *
 *    (0,0) +---------+---------+
 *          |         |         |
 *          |   UL,   |    UR   |
 *          | screen  |         |
 *          +---------+---------+
 *          |         |         |
 *          |   LL    |    LR   |
 *          |         |         |
 *          +---------+---------+ (2 * xsize - 1, 2 * ysize - 1)
 *
 * If panning causes the screen window to move outside the 4 quadrants, swap 
 * quadrants and renormalize the screen coordinates.
 *
 * Example with xsize and ysize = 5. Screen panned to -2, -3 (outside quadrants 
 * in both x and y directions). Swap UL/LR and UR/LL, add xsize to the screen x 
 * coordinate, and add ysize to the screen y coordinate.
 *
 * XXXX = still-valid bitmap data. Blank = new bitmap data that needs to be 
 * calculated.
 *
 * Screen now at 3, 2 after quadrant swap and renormalization.
 * 1. Calculate the update rectangles based on the new screen position
 * 2. Iterate/palette map the update rectangles (XXXX is not recalculated)
 * 3. Blit all 4 quadrant rectangles to the screen (all rectangles now have 
 *    valid bitmap data).
 *
 * (-2, -3)
 *     +---------+
 *     | screen  |
 *     |    +----+----+---------+             +---------+---------+
 *     |    |XXXX|    |         |             |old LR   |old LL   |
 *     +----+----UL   |    UR   |             |    +----+----+    |
 *          |         |         |             |    |    |    |    |
 *          +---------+---------+ --(swap)--> +----+----+----+----+
 *          |         |         |             |    |    |XXXX|    |
 *          |    LL   |    LR   |             |    +----+----+    |
 *          |         |         |             |old UR   |old UL   |
 *          +---------+---------+             +---------+---------+
 *
 * Swappings:
 * if (x < 0 || x > xsize) swap(UL, UR); swap(LL, LR);
 * if (y < 0 || y > ysize) swap(UL, LL); swap(UR, LR);
 *
 * The swap operation swaps only memory pointers and handles, not memory itself.
 *
 * Panning can generate up to 2 update rectangles (1 horizontal, 1 vertical)
 */
static rectangle update_rect[2];

/* The 4 quadrant bitmaps (each of size man_xsize x man_ysize) */
static quadrant quad[4];

#define UL 0	// upper left
#define UR 1	// upper right
#define LL 2	// lower left
#define LR 3	// lower right

/* position of screen (in above coordinate system) */
static int screen_xpos, screen_ypos;

// Constants and variables used in the fast "wave" algorithm

/**
 * Starting values for x and y. Now seems faster to have these as static 
 * globals calculates y = 0, needs dummy line at y = -1 (set to 0's)
 */
static const int wave_ystart[7] = { 3, 1, 3, 1, 0, 1, 0 };
static const int wave_xstart[7] = { 0, 2, 2, 0, 1, 1, 0 };

/* X and Y increments for each wave */
static const int wave_inc[7] = { 4, 4, 4, 4, 2, 2, 2 };

/**
 * Offsets from current pixel location. If all 4 pixels there are equal, set 
 * the current pixel to that value (not used for wave 0). Stored in increasing 
 * pointer order.
 */

static const int wave_xoffs[7][4] = {
	{ 0, 0, 0, 0}, {-2, 2,-2, 2}, {0,-2, 2, 0}, {0,-2, 2, 0},
	{-1, 1,-1, 1}, { 0,-1, 1, 0}, {0,-1, 1, 0} 
};

static const int wave_yoffs[7][4] = {
	{ 0, 0, 0, 0}, {-2,-2, 2, 2}, {-2, 0, 0, 2}, {-2, 0, 0, 2},
	{-1,-1, 1, 1}, {-1, 0, 0, 1}, {-1, 0, 0, 1} 
};

/**
 * Pointer offsets into bitmap, precalculated from above values
 */
static int wave_ptr_offs[7][4];

/**
 * Microsoft syntax for forcing 64-byte alignment (used for aligning 
 * pointstruct arrays in the man_calc_structs).
 */

ALIGNED_(64) man_calc_struct main_man_calc_struct; // used for normal calculation
ALIGNED_(64) man_calc_struct save_man_calc_struct; // used for saving images

/*-------------------------- File/misc functions ----------------------------*/

/**
 * Initialize the settings struct. The name field is what you give in the 
 * logfile to set the setting (not case sensitive). Only the val field can be 
 * modified from logfiles, but default_val can also be set from quickman.cfg. 
 * Here, val should be set the same as the default_val.
 *
 * panrate
 * pan				bitfield; max doesn't matter; autoreset
 * zoomrate
 * zoom				future bitfield; max doesn't matter; autoreset
 * xsize			window functions automatically clip maxes for these; 
 *					autoreset
 * ysize			values less than min size have special meanings
 * maxiters_color	max doesn't really matter, but this only uses 24 bits; 
 *					autoreset.
 * pal_xor			ditto
 * options			bitfield; max doesn't matter
 * spt				bitfield; have to do an external min/max check on this one
 * bst				blit stripe thickness; max doesn't matter
 * pfcmin			10000 * real value
 * pfcmax			10000 * real value
 */
static settings cfg_default = {
	/* name				c_val			d_val			min		max */
	{ "panrate",		DEF_PAN_RATE,	DEF_PAN_RATE,	0,		MAX_PAN_RATE	},
	{ "pan",			0,				0,				0,		0xFFFF			}, 
	{ "zoomrate",		DEF_ZOOM_RATE,	DEF_ZOOM_RATE,	0,		MAX_ZOOM_RATE	},
	{ "zoom",			0,				0,				0,		0xFFFF			}, 
	{ "xsize",			700,			700,			0,		0xFFFF			}, 
	{ "ysize",			700,			700,			0,		0xFFFF			}, 
	{ "maxiters_color", 0,				0,				0,		0xFFFFFF		}, 
	{ "pal_xor",		0,				0,				0,		0xFFFFFF		}, 
	{ "options",		DEF_OPTIONS,	DEF_OPTIONS,	0,		0xFFFF			}, 
	{ "spt",			SPT_DEFAULT,	SPT_DEFAULT,	0,		0xFFFFFF		}, 
	{ "bst",			16,				16,				1,		0xFFFFFF		}, 
	{ "pfcmin",			150,			150,			1,		10000			}, 
	{ "pfcmax",			300,			300,			1,		10000			}, 
};

/**
 * Holds the most recent set of settings read from a file.
 * Copied to the log entry if the entry was valid. Copied to the global config
 * settings if it was read out of quickman.cfg.
 */
static settings cfg_current;

static log_entry* log_entries = NULL;
static int log_pos = 0;
static int log_count = 0;

/**
 * Copy any settings fields that have changed (i.e., are >= 0) from src to dest. 
 * If copy_to_default is 1, also copies changed settings to the default_val 
 * fields (use with quickman.cfg). Then autoreset will reset to the 
 * quickman.cfg default values.
 */
void copy_changed_settings(settings* dest, settings* src, int copy_to_default)
{
	int i;
	setting* s;
	setting* d;

	s = (setting*)src;	// treat structs as arrays
	d = (setting*)dest;
	for (i = 0; i < sizeof(cfg_default) / sizeof(setting); i++)
		if (s[i].c_val >= 0)
		{
			d[i].c_val = s[i].c_val;
			if (copy_to_default)
				d[i].d_val = s[i].c_val;
		}
}

/**
 * Autoreset settings fields to defaults, if so configured. Call before every 
 * image recalculation. Only should be called with global cfg_settings (change 
 * to not take parm?)
 */
void autoreset_settings(settings* dest)
{
	int i;
	setting* d;
	d = (setting*)dest; // treat struct as array
	for (i = 0; i < sizeof(cfg_default) / sizeof(setting); i++)
		if (SETTING_AUTORESET(&d[i]))
			d[i].c_val = d[i].d_val;
}

/* Set all val fields to -1, which invalidates all settings. */
void invalidate_settings(settings* dest)
{
	int i;
	setting* d;
	d = (setting*)dest; // treat struct as array
	for (i = 0; i < sizeof(cfg_default) / sizeof(setting); i++)
		d[i].c_val = -1;
}

/**
 * Read a set of mandelbrot parms into a log entry (if not NULL). This has now 
 * evolved to do a lot more (read in optional config settings and user 
 * commands)
 */
int log_read_entry(log_entry* entry, FILE* fp)
{
	double vals[5];
	int i, j, n, ind, val;
	setting* s, * f;
	unsigned char strs[5][256], * str, c;

	// Initialize cur file settings structure to all invalid (no change)
	invalidate_settings(&cfg_current);

	// Read re, im, mag, iters, pal, and optional commands. To support legacy 
	// logfiles, the five main fields don't need any leading items (they can be 
	// just numbers). Any leading item before a number will be ignored, unless 
	// it's a recognized setting.
	for (i = 0; i < 5;)
	{
		if (feof(fp))
			return 0;
		if (fgets((char*)&strs[i][0], sizeof(strs[0]), fp) == NULL)
			return 0;

		str = &strs[i][0];

		// Skip any leading whitespace
		ind = -1;
		do {
			c = str[++ind];
		} while (c == ' ' || c == '\t'); // null will terminate loop

		// For added robustness, resync on "real", so corrupted files won't get 
		// us out of sync.
		if (!_strnicmp(&str[ind], "real", 4))
			i = 0;

		// Look for any optional commands or settings. This should stay 
		// reasonably fast even with large logfiles.

		s = (setting*)&cfg_default;      // treat structs as arrays
		f = (setting*)&cfg_current;
		for (j = 0; j < sizeof(cfg_default) / sizeof(setting); j++)
			// not case sensitive
			if (!_strnicmp(&str[ind], s[j].name, (size_t)n = strlen(s[j].name)))
			{
				// Can use the function for reading a palette RGB value here 
				// (it will read normal integers too). As a side effect any 24 
				// bit value can be specified as three individual bytes if 
				// desired...
				// Some settings are palette values.

				get_palette_rgb_val(ind + n, str, sizeof(strs[0]), &val);

				// set value if it's within legal range
				if (val >= s[j].min && val <= s[j].max)
					f[j].c_val = val;
				c = 0;  // found a setting; skip the stuff below
				break;
			}
		if (!c)
			continue;

		// Might have an image parameter here (a number with or without leading 
		// items). Strip out any leading non-numeric/non-quote chars, and 
		// ignore comments.
		for (j = ind; j < sizeof(strs[0]); j++)
		{
			if ((c = str[j]) == '/')  // '/' starts a comment
				c = 0;
			if ((c >= '0' && c <= '9') || c == '-' || c == '.' || c == '\"' || !c)
				break;
		}
		if (c)
		{
			// Got something that looks like a number or a " if we get here. 
			// Any bad values will be set to 0.0 (ok). J is long lived (see 
			// below)
			vals[i] = atof(&str[j]);
			i++; // look for next entry
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

		entry->log_settings = cfg_current; // Copy any settings found above

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
			// Get palette from dropdown list. If not found, set to default 
			// palette.
			i = (int)SendDlgItemMessage(hwnd_dialog, 
				IDC_PALETTE, 
				CB_FINDSTRINGEXACT,
				num_builtin_palettes - 1, 
				(LPARAM)&str[j + 1]);

			entry->palette = (i != CB_ERR) ? i : DEFAULT_PAL;
		}
	}
	return 1;
}

/**
 * Open a file for reading. Set bin nonzero to open it in binary mode, else 
 * text mode. 
 * :marked for delete:
 */
FILE* open_file(const char* file, char* msg, int bin)
{
	char s[256];
	
	FILE* fp;
	errno_t err;

	err = fopen_s(&fp, file, bin ? "rb" : "rt");

	if (err != 0)
	{
		if (msg != NULL)
		{
			sprintf_s(s, sizeof(s), "Could not open '%s' for read.%s", 
				file, msg);
			MessageBox(NULL, s, "Warning", 
				MB_OK | MB_ICONWARNING | MB_TASKMODAL);
		}
		return NULL;
	}
	return fp;
}

/**
 * Scan the logfile, dynamically allocate an array for the entries, and fill it 
 * in from the logfile. If init_pos is nonzero, initializes the position to the 
 * beginning.
 */
int log_read(char* file, char* msg, int init_pos)
{
	int i, count;
	FILE* fp;

	log_count = 0;
	if (init_pos)
	{
		log_pos = -1;
		file_tot_time = 0.0; // for benchmarking
	}

	// Kind of inefficient: scan once to get length, then scan again to fill in 
	// array.
	if ((fp = open_file(file, msg, 0)) == NULL)
		return 0;

	for (count = 0; log_read_entry(NULL, fp); count++)
		;

	log_count = count;

	fclose(fp);

	if (!count)
		return 0; // normal cfg files will return here

	if (log_entries != NULL)   // Allocate the array and fill it in
		free(log_entries);
	if ((log_entries = (log_entry*)malloc(count * sizeof(log_entry))) == NULL)
		return 0;

	if ((fp = open_file(file, "", 0)) == NULL)
		return 0;
	for (i = 0; i < count; i++)
		log_read_entry(&log_entries[i], fp);

	fclose(fp);

	return 1;
}

/**
 * Open the logfile for appending and add the current image. Reset position if 
 * reset_pos is 1.
 */
int log_update(char* file, int reset_pos)
{
	char s[512], p[256];
	FILE* fp;
	man_calc_struct* m = &main_man_calc_struct;

	if (fopen_s(&fp, file, "at")) // open for append
	{
		sprintf_s(s, sizeof(s), "Could not open '%s' for write.", file);
		MessageBox(NULL, s, NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
		return 0;
	}

	// For palette, use either number (for builtin palette), or "file" for user 
	// file.
	if (m->palette < num_builtin_palettes)
		sprintf_s(p, sizeof(p), "%d", m->palette);
	else
		sprintf_s(p, sizeof(p), "\"%s\"", palette_file);

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
	return log_read(file, "", reset_pos);
}

/**
 * Get the next or prev entry from the log entry array. Returns the entry.
 */
log_entry* log_get(int next_prevn)
{
	log_entry* e;
	man_calc_struct* m = &main_man_calc_struct;

	if (log_entries == NULL)
		return NULL;

	if (next_prevn)
	{
		if (++log_pos > log_count - 1)
			// log_pos = log_count - 1;	// stop at end
			log_pos = 0;				// wrap to beginning
	}
	else
		if (--log_pos < 0)
			// log_pos = 0;				// stop at beginning
			log_pos = log_count - 1;	// wrap to end

	e = &log_entries[log_pos];

	m->re = e->re;
	m->im = e->im;
	m->mag = e->mag;
	m->max_iters = e->max_iters;
	if (!(status & STAT_PALETTE_LOCKED))
		m->palette = e->palette;

	return e;
}

/**
 * Call this early, before the main window is initialized- could change window 
 * size.
 */
void read_cfg_file(void)
{
	man_calc_struct* m = &main_man_calc_struct;

	// The cfg file is just another logfile, but it shouldn't have any images 
	// in it. If it does, the settings will be reset after each image.
	//
	// This file will update the values in the cfg_settings structure. If it's 
	// missing, the default settings in the structure will be used.

	invalidate_settings(&cfg_current); // initialize all to "no change"

	log_read(CFG_FILE, NULL, 1); // NULL = no error message on failure to read

	// 1 = copy to default_val also
	copy_changed_settings(&cfg_default, &cfg_current, 1);

	// maybe eliminate these separate variables later
	m->xsize = prev_xsize = cfg_default.xsize.c_val;
	m->ysize = prev_ysize = cfg_default.ysize.c_val;

	// just so 1st dialog box doesn't get div by 0
	m->min_dimension = m->xsize > m->ysize ? m->ysize : m->xsize;
	m->max_iters_color = cfg_default.max_iters_color.c_val;
}

/**
 * Find all files ending in .pal and .bmp and add them to the palette dropdown 
 * list. Find all files ending in .log and add to the logfile dropdown list.
 *
 * If calling this more than once (for instance every time the user accesses 
 * the dropdown), should delete old files from the combo box first. 
 */
void add_user_palettes_and_logfiles(void)
{
	HANDLE h = hwnd_dialog;
	int i, n = NUM_ELEM(file_strs);
	LRESULT ind;

	SendDlgItemMessage(h, IDC_PALETTE, CB_DIR, DDL_READONLY | DDL_READWRITE, 
		(LPARAM)"*.pal");
	SendDlgItemMessage(h, IDC_PALETTE, CB_DIR, DDL_READONLY | DDL_READWRITE, 
		(LPARAM)"*.bmp");
	SendDlgItemMessage(h, IDC_LOGFILE, CB_DIR, DDL_READONLY | DDL_READWRITE, 
		(LPARAM)"*.log");

	// Delete any logfiles that were already in the presets list (don't want to 
	// list them twice)
	for (i = 0; i < n; i++)
		if ((ind = SendDlgItemMessage(h, IDC_LOGFILE, CB_FINDSTRINGEXACT, 
			n - 1, (LPARAM)file_strs[i])) >= n)
			SendDlgItemMessage(h, IDC_LOGFILE, CB_DELETESTRING, 
				(WPARAM)ind, 0);
}

/**
 * Timer function: returns current time in TIME_UNITs (changed from previous 
 * versions). Can use two different methods based on #define: 
 * QueryPerformanceCounter is more precise, but has a bug when running on 
 * dual-core CPUs. Randomly selects one or the other core to read the number 
 * from, and gets bogus results if it reads the wrong one.
 */
TIME_UNIT get_timer(void)
{
#ifdef USE_PERFORMANCE_COUNTER
	TIME_UNIT t;
	QueryPerformanceCounter(&t);
	return t;
#else
	return timeGetTime();
#endif
}

/* Get the number of seconds elapsed since start_time. */
double get_seconds_elapsed(TIME_UNIT start_time)
{
	TIME_UNIT t;
#ifdef USE_PERFORMANCE_COUNTER
	QueryPerformanceFrequency(&f);
	QueryPerformanceCounter(&t);
	t.QuadPart -= start_time.QuadPart;
	if (f.QuadPart)
		return (double)t.QuadPart / (double)f.QuadPart;
	return 1e10;
#else
	// Need to use TIME_UNIT (dword) to avoid wrapping issues with timeGetTime().
	// Can subtract in DWORD domain, but not double.
	t = timeGetTime() - start_time;
	return 1e-3 * (double)t;
#endif
}

/**
 * For GetAsyncKeyState: if this bit is set in return value, key is down (MSB 
 * of SHORT)
 */
#define KEYDOWN_BIT		0x8000
#define KEY_LEFT		1
#define KEY_RIGHT		2
#define KEY_UP			4
#define KEY_DOWN		8
#define KEY_CTRL		16
#define KEY_ESC			32
#define KEY_SHIFT		64

/**
 * Return a bitfield indicating the key(s) pressed. Added ASDW as alternate 
 * arrow keys.
 */
int get_keys_pressed(void)
{
	int i, key = 0;
	static SHORT vkeys[] = {
		VK_LEFT,	'A', 
		VK_RIGHT,	'D', 
		VK_UP,		'W', 
		VK_DOWN,	'S', 
		VK_CONTROL, 
		VK_SHIFT 
	};
	static int keybits[] = {
		KEY_LEFT, 
		KEY_LEFT, 
		KEY_RIGHT, 
		KEY_RIGHT, 
		KEY_UP, 
		KEY_UP, 
		KEY_DOWN, 
		KEY_DOWN, 
		KEY_CTRL, 
		KEY_SHIFT 
	};

	for (i = 0; i < NUM_ELEM(vkeys); i++)
		if (GetAsyncKeyState(vkeys[i]) & KEYDOWN_BIT)
			key |= keybits[i];

	return key;
}

/**
 * Reset thread iters accumulators, used to measure thread load balance. Only
 * used for main calculation.
 */
void reset_thread_load_counters(void)
{
	int i;
	man_calc_struct* m = &main_man_calc_struct;

	for (i = 0; i < MAX_THREADS; i++)
		m->thread_states[i].total_iters = 0;
}

/* Reset frames/sec timing values */
void reset_fps_values(void)
{
	total_frames       = 0;
	interval_frames    = 0;
	calc_interval_time = 0.0;
	calc_total_time    = 0.0;
	interval_time      = 0.0;
	total_time         = 0.0;
}

/**
 * Get the real or imaginary coordinate delta based on the pixel delta (offs) 
 * from the image center, the image smaller dimension, and the magnification.
 *
 * Returns value on ST(0), so there should be no precision loss from calling a 
 * function as opposed to doing a macro.
 */
double get_re_im_offs(man_calc_struct* m, long long offs)
{
	return (((double)offs * 4.0) / (double)m->min_dimension) / m->mag;
}

/**
 * Update the image center coordinates (re/im) based on xoffs and yoffs (pixels 
 * from current center). Any time this is done, the pan offsets need to be 
 * reset to 0.

 * Also, call this with pan_xoffs and pan_yoffs to calculate a new re/im from 
 * the current pan offsets (and then reset the offsets).
 */
void update_re_im(man_calc_struct* m, long long xoffs, long long yoffs)
{
	m->re += get_re_im_offs(m, xoffs);
	m->im -= get_re_im_offs(m, yoffs);
	m->pan_xoffs = 0;
	m->pan_yoffs = 0;
}

/**
 * Set the new point and magnification based on x0, x1, y0, y1. If zoom_box 
 * is 0, multiplies/divides the magnification by a fixed value. If zoom_box 
 * is 1, calculates the new zoom from the ratio of the zoom box size (defined 
 * by x0, x1, y0, y1) to the window size. Recenters the image on x0, y0 (if no 
 * zoom box) or the center of the zoom box.

 * Re/im should be updated with the current pan offsets before calling this 
 * (this is done on buttondown events in the window function).
 */
void update_re_im_mag(int zoom_box, int in_outn, int x0, int y0, int x1, int y1)
{
	double tmp_mag, xz, yz;
	int x, y;
	man_calc_struct* m = &main_man_calc_struct;

	tmp_mag = m->mag;

	if (!zoom_box)	// just zoom in or out
	{
		x = x0;     // point center
		y = y0;
		if (in_outn)
			tmp_mag *= MAG_ZOOM_FACTOR;
		else
			tmp_mag *= (1.0 / MAG_ZOOM_FACTOR);
	}
	// Zoom based on box. For zooming in, uses the box size and center to set 
	// the new point and magnification
	else
	{
		x = abs(x0 - x1);
		y = abs(y0 - y1);

		// Get smaller of xzoom and yzoom
		xz = (double)m->xsize / (double)x;
		yz = (double)m->ysize / (double)y;

		if (xz < yz)
			tmp_mag *= xz;
		else
			tmp_mag *= yz;

		x = (x0 + x1) >> 1; // new point center:
		y = (y0 + y1) >> 1; // center of zoom box
	}

	// Get new point center, based on x, y

	//if (in_outn) // uncomment to not center when zooming out
		// re/im already updated with any pan offsets here
		update_re_im(m, x - (m->xsize >> 1), y - (m->ysize >> 1));

	// Update mag
	if (tmp_mag >= MAG_MIN)	// preserve closest min, to allow
		m->mag = tmp_mag;	// zooming back to original mag
}

/* ------------------------ Iteration functions ---------------------------- */

/**
 * Lame unoptimized C iteration function. Just does one point at a time, using 
 * point 0 in the point structure.
 */
static unsigned int iterate_c(man_pointstruct* ps_ptr)
{
	double a, b, x, y, xx, yy, rad;
	unsigned int iters, iter_ct;

	iters = 0;
	iter_ct = ps_ptr->cur_max_iters;

	a = ps_ptr->ab_in[0];
	b = ps_ptr->ab_in[1];
	rad = DIVERGED_THRESH;
	x = y = xx = yy = 0.0;

	do
	{
		y = (x + x) * y + b;
		x = xx - yy + a;
		yy = y * y;
		xx = x * x;
		iters++;
		if ((xx + yy) >= rad)
			break;
	} while (--iter_ct);

	// Store final count and magnitude (squared)

	// use this for tmp storage (mag stored to array below)
	ps_ptr->mag[0] = xx + yy;

	return iters;
}

/* Queue a point to be iterated, for the C iteration function */
static void FASTCALL 
queue_point_c(void* calc_struct, man_pointstruct* ps_ptr, 
	unsigned int* iters_ptr)
{
	man_calc_struct* m;
	unsigned int iters;

	m = (man_calc_struct*)calc_struct;

	iters = m->mandel_iterate(ps_ptr);  // No queuing- just iterate 1 point

	// do this just to match iteration offset of ASM versions
	if (iters != m->max_iters)
		iters++;

	ps_ptr->iters_ptr[0] = iters_ptr;   // Store iters and mag
	*ps_ptr->iters_ptr[0] = iters;
	ps_ptr->iterctr += iters;

	MAG(m, ps_ptr->iters_ptr[0]) = (float)ps_ptr->mag[0];
}

// Code for ASM iteration functions.

// XMM registers for the SSE2 algorithms
#define xmm_x01   xmm0
#define xmm_y01   xmm1
#define xmm_yy01  xmm2
#define xmm_mag   xmm3
#define xmm_x23   xmm4
#define xmm_y23   xmm5
#define xmm_yy23  xmm6
#define xmm_two   xmm7

// XMM registers for the SSE algorithms
#define xmm_x03   xmm0
#define xmm_y03   xmm1
#define xmm_yy03  xmm2
#define xmm_x47   xmm4
#define xmm_y47   xmm5
#define xmm_yy47  xmm6

/**
 * Optimized SSE2 (double precision) ASM algorithm that iterates 4 points at a 
 * time. The point 2,3 calculations are about a half-iteration behind the point 
 * 0,1 calculations, to allow improved interleaving of multiplies and adds. 
 * Additionally, the loop is unrolled twice and divergence detection is done 
 * only every 2nd iteration. After the loop finishes, the code backs out to 
 * check if the points actually diverged during the previous iteration.
 *
 * The divergence detection is done using integer comparisons on the magnitude 
 * exponents (the divergence threshold should be a power of two).
 *
 * 8 total iterations (4 points * 2x unroll) are done per loop. It consists
 * of independent, interleaved iteration blocks:
 *
 * do
 * {
 *    Points 0,1     Points 2,3
 *    y += b;        y *= x;
 *    mag += yy;     x *= x;
 *    x += a;        y *= 2;
 *    yy = y;        mag = x;
 *    yy *= yy;      x -= yy;
 *    y *= x;        y += b;
 *    x *= x;        mag += yy;
 *    y *= 2;        x += a;
 *    mag = x;       yy = y;
 *    x -= yy;       yy *= yy;
 *
 *    [2nd iter: repeat above blocks]
 *
 *    if (any 2nd iter mag diverged)
 *        break;
 * } while (--cur_max_iters);
 *
 * Initial conditions: x = y = yy = 0.0
 *
 * Note: The mags in the first half of the loop are not actually calculated in
 * the loop, but are determined afterwards in the backout calculations from
 * stored intermediate values.
 *
 * From timing measurements, the 53-instruction loop takes ~41 clocks on an 
 * Athlon 64 4000+ 2.4 GHz.
 *
 * Very slow on a Pentium 4 (~100 clocks) due to blocked store fowarding; see 
 * Intel code.
 *
 * sip4
 */
static unsigned int iterate_amd_sse2(man_pointstruct* ps_ptr)
{
	__asm
	{
		// Tried getting rid of the xmm save/restore by having queue_point load 
		// directly into the xmm registers (see qhold_load_xmm.c): overhead 
		// reduction is negligible.

		mov      ebx, ps_ptr		// Get pointstruct pointer. PS4_ macros 
									// below reference [ebx + offset]
		movapd   xmm_x23, PS4_X23	// Restore point states
		movapd   xmm_x01, PS4_X01
		movapd   xmm_two, PS4_TWO
		movapd   xmm_y23, PS4_Y23
		movapd   xmm_y01, PS4_Y01
		movapd   xmm_yy23, PS4_YY23

		mov      edx, DIV_EXP		// Exp for magnitude exponent comparison. 
									// Slower to compare to const directly
		mov      ecx, PS4_CUR_MAX_ITERS // max iters to do this call; always 
										// even
		mov      eax, 0			// iteration counter (for each of the 4 points)
		jmp      skip_top		// Jump past 1st 2 movapds for loop entry, 
								// eliminating
		nop							// need to restore yy01- lower overhead
		nop
		nop							// Achieve the magic alignment (see below)
		nop
		nop
		nop

		// Found that it's important for the end-of-loop branch to be on a 
		// 16-byte boundary (code is slower if not). Choose instructions above 
		// to cause this.
		//
		// v1.0 update: The above no longer holds. Now the magic alignment 
		// seems random...
		// the number of nops above must be determined by trial and error.

	iter_loop:
		movapd   PS4_YY01, xmm_yy01	// save yy01 for mag backout checking
		movapd   PS4_X01, xmm_x01	// save x01 for mag backout checking; 
									// contains xx01 - yy01 here
	skip_top:
		addpd    xmm_y01, PS4_B01	// y01 += b01; faster at top of loop. 
									// Initial y01 = 0
		mulpd    xmm_y23, xmm_x23	/* y23 *= x23; faster at top of loop. */
		add      eax, 2				// update iteration counter; faster here 
									// than 2 insts below
		mulpd    xmm_x23, xmm_x23	/* x23 *= x23 */
		addpd    xmm_x01, PS4_A01	// x01 += a01
		addpd    xmm_y23, xmm_y23	/* y23 *= 2; faster here than mulpd 
									   xmm_y23, xmm_two */
		movapd   xmm_yy01, xmm_y01	// yy01 = y01
		movapd   PS4_X23, xmm_x23	// save xx23 for magnitude backout checking
		mulpd    xmm_yy01, xmm_yy01	/* yy01 *= yy01 */
		subpd    xmm_x23, xmm_yy23	// x23 -= yy23
		mulpd    xmm_y01, xmm_x01	/* y01 *= x01 */
		addpd    xmm_y23, PS4_B23	// y23 += b23
		mulpd    xmm_x01, xmm_x01	/* x01 *= x01 */
		movapd   PS4_YY23, xmm_yy23	// save yy23 for magnitude backout checking
		mulpd    xmm_y01, xmm_two	/* y01 *= 2; add slower here; bb stall */
		addpd    xmm_x23, PS4_A23	// x23 += a23
		movapd   xmm_mag, xmm_x01	// mag01 = x01
		movapd   xmm_yy23, xmm_y23	// yy23 = y23
		subpd    xmm_x01, xmm_yy01	// x01 -= yy01
		mulpd    xmm_yy23, xmm_yy23	/* yy23 *= yy23 */
		addpd    xmm_y01, PS4_B01	// y01 += b01
		mulpd    xmm_y23, xmm_x23	/* y23 *= x23 */
		// ----- Start of 2nd iteration block ------
		addpd    xmm_mag, xmm_yy01
		mulpd    xmm_x23, xmm_x23
		addpd    xmm_x01, PS4_A01
		movapd   xmm_yy01, xmm_y01	// these 2 instrs: faster in this order 
									// than reversed (fixes y23 dep?)
		mulpd    xmm_y23, xmm_two	// (yy01 is apparently just "marked" to get 
									// y01 here; doesn't cause a dep delay)
		movapd   PS4_MAG01, xmm_mag	// save point 0,1 magnitudes for comparison
		movapd   xmm_mag, xmm_x23
		mulpd    xmm_yy01, xmm_yy01
		subpd    xmm_x23, xmm_yy23
		mulpd    xmm_y01, xmm_x01
		addpd    xmm_y23, PS4_B23
		mulpd    xmm_x01, xmm_x01
		addpd    xmm_mag, xmm_yy23
		movapd   PS4_MAG23, xmm_mag	// Save point 2,3 magnitudes. Best here, 
									// despite dep
		cmp      PS4_MEXP0, edx		// Compare the magnitude exponents of 
									// points 0 and 1
		cmovge   ecx, eax			// to the divergence threshold. AMD doesn't 
									// seem to mind
		cmp      PS4_MEXP1, edx	// the store fowarding issue, but Intel does.
		cmovge   ecx, eax			// Conditional moves set ecx to eax on 
									// divergence,
		mulpd    xmm_y01, xmm_two	// breaking the loop.
		addpd    xmm_x23, PS4_A23	// add y, y and mul y, two seem equal speed 
									// here
		movapd   xmm_yy23, xmm_y23
		subpd    xmm_x01, xmm_yy01
		mulpd    xmm_yy23, xmm_yy23
		cmp      PS4_MEXP2, edx		// Compare the magnitude exponents of 
									// points 2 and 3.
		cmovge   ecx, eax
		cmp      PS4_MEXP3, edx
		cmovge   ecx, eax
		cmp      ecx, eax			// Continue iterating until max iters 
									// reached for this call,
		jne      iter_loop			// or one of the points diverged.

		// Exited loop: save iterating point states, and update point iteration 
		// counts. Because divergence detection is only done every 2 iterations, 
		// need to "back out" and see if a point diverged the previous iteration. 
		// Calculate previous magnitudes from stored values and save in expp. 
		// Caller can then use DIVERGED and DIVERGED_PREV macros to detect 
		// if/when the point diverged.

		// Structure contents here: 
		//		PS4.x = xx01 - yy01;
		//		PS4.yy = yy01; 
		//		PS4.x + 16 = xx23; 
		//		PS4.yy + 16 = yy23
		// Order here seems fastest, despite dependencies
		// Really only need to calculate prev mags for points that diverged... 
		// could save overhead

		movapd   PS4_Y01, xmm_y01	// save y01 state
		movapd   PS4_Y23, xmm_y23	// save y23 state

		mulpd    xmm_two, PS4_YY01	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy01 */
		movapd   xmm_mag, PS4_X23	// tmp2 = xx23
		addpd    xmm_two, PS4_X01	/* get mag01 = xx01 - yy01 + 2 * yy01 = 
									   xx01 + yy01 */
		addpd    xmm_mag, PS4_YY23	// get mag23 = xx23 + yy23
		movapd   PS4_MAGPREV01, xmm_two		// store prev_mag 01
		movapd   PS4_MAGPREV23, xmm_mag		// store prev_mag 23

		xor ecx, ecx				// Get a 0
		add      PS4_ITERCTR_L, eax	// Update iteration counter. Multiply by 36 
									// to get effective flops.
		adc      PS4_ITERCTR_H, ecx	// Update iterctr high dword

		movapd   PS4_YY01, xmm_yy01	// save yy01 state
		movapd   PS4_X23, xmm_x23	// save x23 state
		movapd   PS4_X01, xmm_x01	// save x01 state
		movapd   PS4_YY23, xmm_yy23	// save yy23 state

		add      PS4_ITERS0, eax	// update point iteration counts
		add      PS4_ITERS1, eax
		add      PS4_ITERS2, eax
		add      PS4_ITERS3, eax

		// return value (iters done per point) is in eax

		// 0xF3 prefix for AMD single byte return fix: does seem to make a 
		// difference
		// Don't put any code between here and ret. Watch for compiler-inserted 
		// pops if using extra registers.

		//pop      ebx				// compiler pushes ebx, ebp on entry
		//pop      ebp				// accessing ebp gives compiler warning
		//__emit(0xF3);
		//ret
	}
}

/**
 * Intel version. Basically the same as the AMD version, but doesn't do 
 * magnitude comparison in the INT domain because this causes a blocked 
 * store-forwarding penalty (from XMM store to INT load) of 40-50 clock cycles. 
 * Uses cmppd/movmskpd instead.
 *
 * Changed to use same initial conditions as AMD code. No consistent speed 
 * difference vs. old version. Appears to be about 0.3% slower on the home 
 * benchmark, but 0.5% faster on bmark.log.

 * From timing measurements, the loop takes 50 clocks on a Pentium D 820 2.8 
 * GHz.
 *
 * sip4
 */
static unsigned int iterate_intel_sse2(man_pointstruct* ps_ptr)
{
	__asm
	{
		mov      ebx, ps_ptr		// Get pointstruct pointer
		movapd   xmm_yy01, PS4_YY01	// Restore point states
		movapd   xmm_y01, PS4_Y01
		movapd   xmm_x23, PS4_X23
		movapd   xmm_x01, PS4_X01
		movapd   xmm_two, PS4_TWO
		movapd   xmm_y23, PS4_Y23
		movapd   xmm_yy23, PS4_YY23
		addpd    xmm_y01, PS4_B01	// pre-add y01 to get correct initial 
									// condition

		mov      eax, 2			// iteration counter (for each of the 4 points)
		jmp      skip_top

	iter_loop:						// alignment doesn't seem to matter on Intel
		movapd   PS4_YY01, xmm_yy01	// save yy01 for mag backout checking
		movapd   PS4_X01, xmm_x01	// save x01 for mag backout checking; 
									// contains xx01 - yy01 here
		add      eax, 2				// update iteration counter
	skip_top:
		mulpd    xmm_x23, xmm_x23	/* x23 *= x23 */
		addpd    xmm_x01, PS4_A01	// x01 += a01
		addpd    xmm_y23, xmm_y23	/* y23 *= 2; faster here than mulpd 
									   xmm_y23, xmm_two */
		movapd   xmm_yy01, xmm_y01	// yy01 = y01
		movapd   PS4_X23, xmm_x23	// save xx23 for magnitude backout checking
		mulpd    xmm_yy01, xmm_yy01	/* yy01 *= yy01 */
		subpd    xmm_x23, xmm_yy23	// x23 -= yy23
		mulpd    xmm_y01, xmm_x01	/* y01 *= x01 */
		addpd    xmm_y23, PS4_B23	// y23 += b23
		mulpd    xmm_x01, xmm_x01	/* x01 *= x01 */
		movapd   PS4_YY23, xmm_yy23	// save yy23 for magnitude backout checking
		mulpd    xmm_y01, xmm_two	/* y01 *= 2; add slower here */
		addpd    xmm_x23, PS4_A23	// x23 += a23
		movapd   xmm_mag, xmm_x01	// mag01 = x01
		movapd   xmm_yy23, xmm_y23	// yy23 = y23
		subpd    xmm_x01, xmm_yy01	// x01 -= yy01
		mulpd    xmm_yy23, xmm_yy23	/* yy23 *= yy23 */
		addpd    xmm_y01, PS4_B01	// y01 += b01
		mulpd    xmm_y23, xmm_x23	/* y23 *= x23 */
		// ----- Start of 2nd iteration block ------
		addpd    xmm_mag, xmm_yy01
		mulpd    xmm_x23, xmm_x23
		addpd    xmm_x01, PS4_A01
		mulpd    xmm_y23, xmm_two
		movapd   xmm_yy01, xmm_y01
		movapd   PS4_MAG01, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltpd xmm_mag, PS4_RAD	// compare point 0, 1 magnitudes 
									// (mag >= rad): let cpu reorder these
		movmskpd edx, xmm_mag		// save result in edx
		movapd   xmm_mag, xmm_x23
		mulpd    xmm_yy01, xmm_yy01
		subpd    xmm_x23, xmm_yy23
		mulpd    xmm_y01, xmm_x01
		addpd    xmm_y23, PS4_B23
		mulpd    xmm_x01, xmm_x01
		addpd    xmm_mag, xmm_yy23
		mulpd    xmm_y01, xmm_two
		addpd    xmm_x23, PS4_A23
		add      edx, edx		// shift point 01 mag compare results left 2
		add      edx, edx
		movapd   xmm_yy23, xmm_y23
		subpd    xmm_x01, xmm_yy01
		movapd   PS4_MAG23, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltpd xmm_mag, PS4_RAD	// compare point 2, 3 magnitudes
		mulpd    xmm_yy23, xmm_yy23
		addpd    xmm_y01, PS4_B01
		movmskpd ecx, xmm_mag
		mulpd    xmm_y23, xmm_x23
		or ecx,	edx					// Continue iterating until max iters 
									// reached for this call,
		jnz      done				// or one of the points diverged.
		cmp      PS4_CUR_MAX_ITERS, eax		// No penalty for comparing from 
											// memory vs. register here
		jne      iter_loop

	done:
		subpd    xmm_y01, PS4_B01	// subtract out pre-add (see loop top)
		movapd   PS4_Y01, xmm_y01	// save y01 state
		movapd   PS4_Y23, xmm_y23	// save y23 state

		// Get previous magnitudes. See AMD code
		mulpd    xmm_two, PS4_YY01	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy01 */
		movapd   xmm_mag, PS4_X23	// tmp2 = xx23
		addpd    xmm_two, PS4_X01	/* get mag01 = xx01 - yy01 + 2 * yy01 = 
									   xx01 + yy01 */
		addpd    xmm_mag, PS4_YY23	// get mag23 = xx23 + yy23
		movapd   PS4_MAGPREV01, xmm_two // store prev_mag 01
		movapd   PS4_MAGPREV23, xmm_mag // store prev_mag 23

		xor edx, edx				// Get a 0
		add      PS4_ITERCTR_L, eax	// Update iteration counter. Multiply by 36 
									// to get effective flops.
		adc      PS4_ITERCTR_H, edx	// Update iterctr high dword

		movapd   PS4_YY01, xmm_yy01	// save yy01 state
		movapd   PS4_X23, xmm_x23	// save x23 state
		movapd   PS4_X01, xmm_x01	// save x01 state
		movapd   PS4_YY23, xmm_yy23	// save yy23 state

		add      PS4_ITERS0, eax	// update point iteration counts
		add      PS4_ITERS1, eax
		add      PS4_ITERS2, eax
		add      PS4_ITERS3, eax
	}
}

/**
 * SSE (single precision) algorithm for AMD; iterates 8 points at a time,
 * or 16 iterations per loop. Based on AMD SSE2 algorithm.
 *
 * Loop appears to take 42.25 clocks- the 8 extra int compare/cmov
 * instructions add 2 extra clocks per loop vs. the SSE2 algorithm.
 * Tried rearranging them but current order seems best.
 *
 * sip8
 */
static unsigned int iterate_amd_sse(man_pointstruct* ps_ptr)
{
	__asm
	{
		mov      ebx, ps_ptr		// Get pointstruct pointer
		movaps   xmm_x47, PS8_X47	// Restore point states
		movaps   xmm_x03, PS8_X03
		movaps   xmm_two, PS8_TWO
		movaps   xmm_y47, PS8_Y47
		movaps   xmm_y03, PS8_Y03
		movaps   xmm_yy47, PS8_YY47

		mov      edx, DIV_EXP_FLOAT	// Exp for magnitude exponent comparison. 
									// Slower to compare to const directly
		mov      ecx, PS8_CUR_MAX_ITERS	// max iters to do this call; always 
										// even
		mov      eax, 0		// Iteration counter (for each of the 4 points)
		jmp      skip_top			// Allows removing yy03 restore above
		nop
		nop							// Achieve the magic alignment...
		nop
		nop
		nop
		nop
		nop
		nop
		nop

	iter_loop:
		movaps   PS8_YY03, xmm_yy03	// save yy03 for mag backout checking
		movaps   PS8_X03, xmm_x03	// save x03 for mag backout checking; 
									// contains xx03 - yy03 here
	skip_top:
		addps    xmm_y03, PS8_B03	// y03 += b03; faster at top of loop. 
									// Initial y03 = 0
		mulps    xmm_y47, xmm_x47	/* y47 *= x47; faster at top of loop. */
		add      eax, 2				// update iteration counter; faster here 
									// than 2 insts below
		mulps    xmm_x47, xmm_x47	// x47 *= x47
		addps    xmm_x03, PS8_A03	// x03 += a03
		addps    xmm_y47, xmm_y47	/* y47 *= 2; faster here than mulps 
									   xmm_y47, xmm_two */
		movaps   xmm_yy03, xmm_y03	// yy03 = y03
		movaps   PS8_X47, xmm_x47	// save xx47 for magnitude backout checking
		mulps    xmm_yy03, xmm_yy03	/* yy03 *= yy03 */
		subps    xmm_x47, xmm_yy47	// x47 -= yy47
		mulps    xmm_y03, xmm_x03	/* y03 *= x03 */
		addps    xmm_y47, PS8_B47	// y47 += b47
		mulps    xmm_x03, xmm_x03	/* x03 *= x03 */
		movaps   PS8_YY47, xmm_yy47	// save yy47 for magnitude backout checking
		mulps    xmm_y03, xmm_two	/* y03 *= 2; add slower here; bb stall */
		addps    xmm_x47, PS8_A47	// x47 += a47
		movaps   xmm_mag, xmm_x03	// mag03 = x03
		movaps   xmm_yy47, xmm_y47	// yy47 = y47
		subps    xmm_x03, xmm_yy03	// x03 -= yy03
		mulps    xmm_yy47, xmm_yy47	/* yy47 *= yy47 */
		addps    xmm_y03, PS8_B03	// y03 += b03
		mulps    xmm_y47, xmm_x47	/* y47 *= x47 */
		// ----- Start of 2nd iteration block ------
		addps    xmm_mag, xmm_yy03
		mulps    xmm_x47, xmm_x47
		addps    xmm_x03, PS8_A03
		movaps   xmm_yy03, xmm_y03	// these 2 instrs: faster in this order 
									// than reversed (fixes y47 dep?)
		mulps    xmm_y47, xmm_two	// (yy03 is apparently just "marked" to get 
									// y03 here; doesn't cause a dep delay)
		movaps   PS8_MAG03, xmm_mag	// save point 0-3 magnitudes for comparison
		movaps   xmm_mag, xmm_x47
		mulps    xmm_yy03, xmm_yy03
		subps    xmm_x47, xmm_yy47
		mulps    xmm_y03, xmm_x03
		addps    xmm_y47, PS8_B47
		mulps    xmm_x03, xmm_x03
		addps    xmm_mag, xmm_yy47
		movaps   PS8_MAG47, xmm_mag	// Save point 4-7 magnitudes. Best here, 
									// despite dep
		cmp      PS8_MEXP0, edx		// Compare the magnitude exponents of 
									// points 0-3
		cmovge   ecx, eax			// to the divergence threshold. AMD doesn't 
									// seem to mind
		cmp      PS8_MEXP1, edx		// the store fowarding issue, but Intel 
									// does.
		cmovge   ecx, eax			// Conditional moves set ecx to eax on 
									// divergence,
		cmp      PS8_MEXP2, edx		// breaking the loop.
		cmovge   ecx, eax
		cmp      PS8_MEXP3, edx
		cmovge   ecx, eax
		mulps    xmm_y03, xmm_two	// add y, y and mul y, two seem equal speed 
									// here
		addps    xmm_x47, PS8_A47
		movaps   xmm_yy47, xmm_y47
		subps    xmm_x03, xmm_yy03
		mulps    xmm_yy47, xmm_yy47
		cmp      PS8_MEXP4, edx		// Compare the magnitude exponents of 
									// points 4-7.
		cmovge   ecx, eax
		cmp      PS8_MEXP5, edx
		cmovge   ecx, eax
		cmp      PS8_MEXP6, edx
		cmovge   ecx, eax
		cmp      PS8_MEXP7, edx
		//cmovge   ecx, eax			// jge done seems a hair faster. When 
									// changing instrs, don't forget
		jge      done				// to adjust nops to put jne iter_loop on a 
									// 16-byte boundary.
		cmp      ecx, eax			// Continue iterating until max iters 
									// reached for this call,
		jne      iter_loop			// or one of the points diverged.

	done :
		// Get previous magnitudes. See AMD SSE2 code
		movaps   PS8_Y03, xmm_y03	// save y03 state
		movaps   PS8_Y47, xmm_y47	// save y47 state

		mulps    xmm_two, PS8_YY03	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy03 */
		movaps   xmm_mag, PS8_X47	// tmp2 = xx47
		addps    xmm_two, PS8_X03	/* get mag03 = xx03 - yy03 + 2 * yy03 = 
									   xx03 + yy03 */
		addps    xmm_mag, PS8_YY47	// get mag47 = xx47 + yy47
		movaps   PS8_MAGPREV03, xmm_two		// store prev_mag 03
		movaps   PS8_MAGPREV47, xmm_mag		// store prev_mag 47

		xor ecx, ecx				// Get a 0
		add      PS8_ITERCTR_L, eax	// Update iteration counter. Multiply by 72 
									// to get effective flops.
		adc      PS8_ITERCTR_H, ecx	// Update iterctr high dword

		movaps   PS8_YY03, xmm_yy03	// save yy03 state
		movaps   PS8_X47, xmm_x47	// save x47 state
		movaps   PS8_X03, xmm_x03	// save x03 state
		movaps   PS8_YY47, xmm_yy47	// save yy47 state

		add      PS8_ITERS0, eax	// update point iteration counts
		add      PS8_ITERS1, eax
		add      PS8_ITERS2, eax
		add      PS8_ITERS3, eax
		add      PS8_ITERS4, eax
		add      PS8_ITERS5, eax
		add      PS8_ITERS6, eax
		add      PS8_ITERS7, eax
		// return value (iterations done per point) is in eax
	}
}

/**
 * SSE (single precision) algorithm for Intel; iterates 8 points at a time, or 
 * 16 iterations per loop. Based on Intel SSE2 algorithm.
 *
 * sip8
 */
static unsigned int iterate_intel_sse(man_pointstruct* ps_ptr)
{
	__asm
	{
		mov      ebx, ps_ptr		// Get pointstruct pointer
		movaps   xmm_yy03, PS8_YY03	// Restore point states
		movaps   xmm_x47, PS8_X47
		movaps   xmm_x03, PS8_X03
		movaps   xmm_two, PS8_TWO
		movaps   xmm_y47, PS8_Y47
		movaps   xmm_y03, PS8_Y03
		movaps   xmm_yy47, PS8_YY47
		addps    xmm_y03, PS8_B03	// pre-add y03 to get correct initial 
									//condition
		mov      eax, 2			// iteration counter (for each of the 4 points)
		jmp      skip_top

	iter_loop:					// alignment doesn't seem to matter on Intel
		movaps   PS8_YY03, xmm_yy03	// save yy03 for mag backout checking
		movaps   PS8_X03, xmm_x03	// save x03 for mag backout checking; 
									// contains xx03 - yy03 here
		add      eax, 2				// update iteration counter
	skip_top :
		mulps    xmm_x47, xmm_x47	/* x47 *= x47 */
		addps    xmm_x03, PS8_A03	// x03 += a03
		addps    xmm_y47, xmm_y47	/* y47 *= 2; faster here than mulps 
									   xmm_y47, xmm_two */
		movaps   xmm_yy03, xmm_y03	// yy03 = y03
		movaps   PS8_X47, xmm_x47	// save xx47 for magnitude backout checking
		mulps    xmm_yy03, xmm_yy03	/* yy03 *= yy03 */
		subps    xmm_x47, xmm_yy47	// x47 -= yy47
		mulps    xmm_y03, xmm_x03	/* y03 *= x03 */
		addps    xmm_y47, PS8_B47	// y47 += b47
		mulps    xmm_x03, xmm_x03	/* x03 *= x03 */
		movaps   PS8_YY47, xmm_yy47	// save yy47 for magnitude backout checking
		mulps    xmm_y03, xmm_two	/* y03 *= 2; add slower here */
		addps    xmm_x47, PS8_A47	// x47 += a47
		movaps   xmm_mag, xmm_x03	// mag03 = x03
		movaps   xmm_yy47, xmm_y47	// yy47 = y47
		subps    xmm_x03, xmm_yy03	// x03 -= yy03
		mulps    xmm_yy47, xmm_yy47	/* yy47 *= yy47 */
		addps    xmm_y03, PS8_B03	// y03 += b03
		mulps    xmm_y47, xmm_x47	/* y47 *= x47 */
		// ----- Start of 2nd iteration block ------
		addps    xmm_mag, xmm_yy03
		mulps    xmm_x47, xmm_x47
		addps    xmm_x03, PS8_A03
		mulps    xmm_y47, xmm_two
		movaps   xmm_yy03, xmm_y03
		
		movapd   PS8_MAG03, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltps xmm_mag, PS8_RAD	// compare point 0-3 magnitudes 
									// (mag >= rad): let cpu reorder these
		movmskps edx, xmm_mag		// save result in edx
		movaps   xmm_mag, xmm_x47
		mulps    xmm_yy03, xmm_yy03
		subps    xmm_x47, xmm_yy47
		mulps    xmm_y03, xmm_x03
		addps    xmm_y47, PS8_B47
		mulps    xmm_x03, xmm_x03
		addps    xmm_mag, xmm_yy47
		mulps    xmm_y03, xmm_two
		addps    xmm_x47, PS8_A47
		shl      edx, 4			// shift point 0-3 mag compare results left 4
		movaps   xmm_yy47, xmm_y47
		subps    xmm_x03, xmm_yy03
		movapd   PS8_MAG47, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltps xmm_mag, PS8_RAD	// compare point 4-7 magnitudes
		mulps    xmm_yy47, xmm_yy47
		addps    xmm_y03, PS8_B03
		movmskps ecx, xmm_mag
		mulps    xmm_y47, xmm_x47
		or ecx, edx		// Continue iterating until max iters reached for this 
						// call,
		jnz      done	// or one of the points diverged.
		cmp      PS8_CUR_MAX_ITERS, eax		// No penalty for comparing from 
											// memory vs. register here
		jne      iter_loop

	done :
		subps    xmm_y03, PS8_B03	// subtract out pre-add (see loop top)
		movaps   PS8_Y03, xmm_y03	// save y03 state
		movaps   PS8_Y47, xmm_y47	// save y47 state

		// Get previous magnitudes. See AMD SSE2 code
		mulps    xmm_two, PS8_YY03	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy03 */
		movaps   xmm_mag, PS8_X47	// tmp2 = xx47
		addps    xmm_two, PS8_X03	/* get mag03 = xx03 - yy03 + 2 * yy03 = 
									   xx03 + yy03 */
		addps    xmm_mag, PS8_YY47	// get mag47 = xx47 + yy47
		movaps   PS8_MAGPREV03, xmm_two	// store prev_mag 03
		movaps   PS8_MAGPREV47, xmm_mag	// store prev_mag 47

		xor edx, edx				// Get a 0
		add      PS8_ITERCTR_L, eax	// Update iteration counter. Multiply 
									// by 72 to get effective flops.
		adc      PS8_ITERCTR_H, edx	// Update iterctr high dword

		movaps   PS8_YY03, xmm_yy03	// save yy03 state
		movaps   PS8_X47, xmm_x47	// save x47 state
		movaps   PS8_X03, xmm_x03	// save x03 state
		movaps   PS8_YY47, xmm_yy47	// save yy47 state

		add      PS8_ITERS0, eax	// update point iteration counts
		add      PS8_ITERS1, eax
		add      PS8_ITERS2, eax
		add      PS8_ITERS3, eax
		add      PS8_ITERS4, eax
		add      PS8_ITERS5, eax
		add      PS8_ITERS6, eax
		add      PS8_ITERS7, eax
	}
}

/* --------------------------- Queuing functions --------------------------- */

/**
 * The queue_status field of the pointstruct structure keeps track of which 
 * point queue slots are free. Each 3 bits gives a free slot.
 * Push: shift left 3, OR in free slot number.
 * Pop: slot number = low 3 bits, shift right 3.
 * Queue is full when queue_status == QUEUE_FULL.
 * Initialize with 3, 2, 1, 0 in bits 11-0, QUEUE_FULL in bits 15-12
 * Also used for single precision algorithm, to track up to 8 queue slots
 */
#define QUEUE_FULL   0xF

// Condition indicating point[ind] diverged
#define DIVERGED(p, ind)		(((int *)p->mag)[1 + (ind << 1)] >= DIV_EXP)

// Condition indicating point[ind] diverged on the previous iteration
#define DIVERGED_PREV(p, ind)	(((int *)p->magprev)[1 + (ind << 1)] >= DIV_EXP)

// For single-precision (SSE) version
#define DIVERGED_S(p, ind)		(((int *)p->mag)[ind] >= DIV_EXP_FLOAT)
#define DIVERGED_PREV_S(p, ind)	(((int *)p->magprev)[ind] >= DIV_EXP_FLOAT)

/**
 * Queue a point for iteration using the 4-point SSE2 algorithm. On entry, 
 * ps_ptr->ab_in should contain the real and imaginary parts of the point to 
 * iterate on, and iters_ptr should be the address of the point in the 
 * iteration count array.
 *
 * Can't assume all 4 point b's will be the same (pixels on the previous line 
 * may still be iterating).
 *
 * Rewriting this in ASM would probably save a lot of overhead (important for
 * realtime zooming and panning)
 *
 * sqp
 */
static void FASTCALL 
queue_4point_sse2(void* calc_struct, man_pointstruct* ps_ptr, 
	unsigned int* iters_ptr)
{
	unsigned int i, iters, max, queue_status, * ptr;
	man_calc_struct* m;

	m = (man_calc_struct*)calc_struct;

	queue_status = ps_ptr->queue_status;

	// If all points in use, iterate to clear at least one point first
	if (queue_status == QUEUE_FULL)
	{
		// Returns (iters done) if any point hit max iterations, or diverged
		m->mandel_iterate(ps_ptr);

		max = 0;
		for (i = 0; i < 4; i++) // compiler fully unrolls this
		{
			// Find which point(s) are done and retire them (store iteration 
			// count to array). Iteration counts will later be mapped to colors 
			// using the current palette. Timing test: removing array stores 
			// results in NO time savings on bmark.log.

			iters = ps_ptr->iters[i];
			if (DIVERGED(ps_ptr, i))
			{
				// If actually diverged on previous iteration, dec iters.
				// *ps_ptr->iters_ptr[i] = iters - DIVERGED_PREV(ps_ptr, i);

				// This is about 3% slower (with the branch and extra stores) 
				// than the old code (above).
				// Needed to get the correct magnitude for the normalized 
				// iteration count algorithm.

				ptr = ps_ptr->iters_ptr[i];
				if (DIVERGED_PREV(ps_ptr, i))
				{
					*ptr = iters - 1;
					MAG(m, ptr) = (float)ps_ptr->magprev[i];
				}
				else
				{
					*ptr = iters;
					MAG(m, ptr) = (float)ps_ptr->mag[i];
				}

				// Push free slot. Use this form to allow compiler to use the 
				// lea instruction.
				queue_status = queue_status * 8 + i;
			}
			// Gets here most often. See if this point has the most accumulated 
			// iterations. Also check if point reached max iters and retire if 
			// so. Definite overhead improvement to combine the max iters check 
			// with the max check- measurable with small max_iters (e.g., when 
			// realtime zooming).
			else
				if (iters >= max)
					if (iters == m->max_iters)
					{
						// don't need mag store for max_iters
						*ps_ptr->iters_ptr[i] = iters;

						queue_status = queue_status * 8 + i; // Push free slot
					}
					else
						max = iters;

		}
		// Set the maximum iterations to do next loop: max iters - iters 
		// already done. The next loop must break if the point with the most 
		// accumulated iterations (max) reaches max_iters.
		ps_ptr->cur_max_iters = m->max_iters - max;

		// Most common case: comes here with one point free. Retiring multiple 
		// points only happens about 1-5% of the time for complex images. For 
		// images with vast areas of a single color, can go to 50%.
	}

	i = queue_status & 3;                     // Get next free slot
	ps_ptr->queue_status = queue_status >> 3; // Pop free slot

	// Initialize pointstruct fields
	ps_ptr->a[i] = ps_ptr->ab_in[0];	// Set input point
	ps_ptr->b[i] = ps_ptr->ab_in[1];
	ps_ptr->y[i] = 0.0;					// Set initial conditions
	ps_ptr->x[i] = 0.0;
	ps_ptr->yy[i] = 0.0;
	ps_ptr->iters[i] = 0;
	ps_ptr->iters_ptr[i] = iters_ptr;
}

/**
 * Similar queuing function for the 8-point SSE algorithm.
 *
 * sqp8
 */
static void FASTCALL 
queue_8point_sse(void* calc_struct, man_pointstruct* ps_ptr, 
	unsigned* iters_ptr)
{
	unsigned i, iters, max, queue_status, * ptr;
	man_calc_struct* m;

	m = (man_calc_struct*)calc_struct;
	queue_status = ps_ptr->queue_status;

	if (queue_status == QUEUE_FULL)
	{
		m->mandel_iterate(ps_ptr);

		max = 0;
		for (i = 0; i < 8; i++)
		{
			iters = ps_ptr->iters[i];
			if (DIVERGED_S(ps_ptr, i))
			{
				/**ps_ptr->iters_ptr[i] = iters - DIVERGED_PREV_S(ps_ptr, i);*/

				ptr = ps_ptr->iters_ptr[i];
				if (DIVERGED_PREV_S(ps_ptr, i))
				{
					*ptr = iters - 1;
					MAG(m, ptr) = ((float*)ps_ptr->magprev)[i];
				}
				else
				{
					*ptr = iters;
					MAG(m, ptr) = ((float*)ps_ptr->mag)[i];
				}
				queue_status = queue_status * 8 + i;
			}
			else
				if (iters >= max)
					if (iters == m->max_iters)
					{
						*ps_ptr->iters_ptr[i] = iters;
						queue_status = queue_status * 8 + i;
					}
					else
						max = iters;
		}
		ps_ptr->cur_max_iters = m->max_iters - max;
	}

	i = queue_status & 7;
	ps_ptr->queue_status = queue_status >> 3;

	// Initialize pointstruct fields as packed 32-bit floats
	((float*)ps_ptr->a)[i] = (float)ps_ptr->ab_in[0];	// Set input point- 
														// convert from doubles
	((float*)ps_ptr->b)[i] = (float)ps_ptr->ab_in[1];	// generated by the 
														// main loop
	((float*)ps_ptr->y)[i] = 0.0;	// Set initial conditions
	((float*)ps_ptr->x)[i] = 0.0;
	((float*)ps_ptr->yy)[i] = 0.0;
	ps_ptr->iters[i] = 0;
	ps_ptr->iters_ptr[i] = iters_ptr;
}

/**
 * Calculate the image, using the currently set precision and algorithm. 
 * Calculations in here are always done in double (or extended) precision, 
 * regardless of the iteration algorithm's precision.
 *
 * Now called from multiple threads. Calculates a list of stripes from the 
 * thread state structure (passed in PARAM). See man_calculate().
 *
 * smc
 */
unsigned int CALLBACK man_calculate_threaded(LPVOID param)
{
	int i, n, x, y, xstart, xend, ystart, yend, line_size, points_guessed;
	unsigned int* iters_ptr;
	man_pointstruct* ps_ptr;
	thread_state* t;
	stripe* s;
	man_calc_struct* m;

	t = (thread_state*)param;
	s = t->stripes;
	n = t->num_stripes;
	ps_ptr = t->ps_ptr;
	m = (man_calc_struct*)t->calc_struct;

	line_size = m->iter_data_line_size;
	points_guessed = 0;

	// Calculate all the stripes. Needs to handle num_stripes == 0
	for (i = 0; i < n; i++)
	{
		// Use these for benchmarking thread creation/execution overhead
		// return 0;                  // for CreateThread method
		// SetEvent(t->done_event);   // for QueueUserWorkItem method
		// return 0;

		xstart = s->xstart;
		xend = s->xend;
		ystart = s->ystart;
		yend = s->yend;

		// Optimization for panning: set alg to exact mode for very thin 
		// regions. Due to the Fast algorithm's 4x4 cell size it often computes 
		// more pixels than Exact for these regions. Effect is most apparent 
		// with high iter count images.

#define FE_SWITCHOVER_THRESH  2 // only do it for 1-pixel wide regions for now. 
								// best value TBD...

		m->cur_alg = m->alg;
		if ((xend - xstart) < FE_SWITCHOVER_THRESH || 
			(yend - ystart) < FE_SWITCHOVER_THRESH)
			m->cur_alg |= ALG_EXACT;

		// Main loop. Queue each point in the image for iteration. Queue_point 
		// will return immediately if its queue isn't full (needs 4 points for 
		// the asm version), otherwise it will iterate on all the points in the 
		// queue.

		if (m->cur_alg & ALG_EXACT) // Exact algorithm: calculates every pixel
		{
			y = ystart;
			do
			{
				x = xstart;

				// Load IM coordinate from the array
				ps_ptr->ab_in[1] = m->img_im[y];
				
				iters_ptr = m->iter_data + y * line_size + x;
				do
				{
					// Load RE coordinate from the array
					ps_ptr->ab_in[0] = m->img_re[x];
					
					m->queue_point(m, ps_ptr, iters_ptr++);
				} while (++x <= xend);
			} while (++y <= yend);
		}
		else // Fast "wave" algorithm from old code: guesses pixels.
		{
			int wave, xoffs, inc, p0, p1, p2, p3, offs0, offs1, offs2, offs3;

			// Doing the full calculation (all waves) on horizontal chunks to 
			// improve cache locality gives no speedup (tested before realtime 
			// zooming was implemented- maybe should test again).

			for (wave = 0; wave < 7; wave++)
			{
				inc = wave_inc[wave];
				y = wave_ystart[wave] + ystart;

				// Special case for wave 0 (always calculates all pixels). 
				// Makes realtime zooming measurably faster. X starts at xstart 
				// for wave 0, so can use do-while. For Y, need to calculate 
				// all waves even if out of range, because subsequent waves 
				// look forward to pixels calculated in previous waves (wave 0 
				// starts at y = 3).

				// it's faster with the special case inside the wave loop than 
				// outside
				if (!wave)
				{
					do
					{
						x = xstart;

						// Load IM coordinate from the array
						ps_ptr->ab_in[1] = m->img_im[y];
						
						// adding a line to the ptr every y loop is slower
						iters_ptr = m->iter_data + y * line_size + x;

						do {
							// Load RE coordinate from the array
							ps_ptr->ab_in[0] = m->img_re[x];
							m->queue_point(m, ps_ptr, iters_ptr);
							iters_ptr += inc;
							x += inc;
						} while (x <= xend);
					} while ((y += inc) <= yend);
				}
				else  // waves 1-6 check neighboring pixels
				{
					// pointer offsets of neighboring pixels
					offs0 = wave_ptr_offs[wave][0];
					offs1 = wave_ptr_offs[wave][1];
					offs2 = wave_ptr_offs[wave][2];
					offs3 = wave_ptr_offs[wave][3];

					xoffs = wave_xstart[wave] + xstart;

					do
					{
						x = xoffs;
						ps_ptr->ab_in[1] = m->img_im[y];
						iters_ptr = m->iter_data + y * line_size + x;

						// No faster to have a special case for waves 1 and 4 
						// that loads only 2 pixels/loop
						while (x <= xend)
						{
							// If all 4 neighboring pixels (p0 - p3) are the 
							// same, set this pixel to their value, else 
							// iterate.

							p0 = iters_ptr[offs0];
							p1 = iters_ptr[offs1];
							p2 = iters_ptr[offs2];
							p3 = iters_ptr[offs3];

							// can't use sum compares here (causes corrupted 
							// pixels)
							if (p0 == p1 && p0 == p2 && p0 == p3)
							{
								// aargh... compiler (or AMD CPU) generates 
								// different performance on zoomtest depending 
								// on which point is stored here. They're all 
								// the same...
								// p3: 18.5s  p2: 18.3s  p1: 18.7s  p0: 19.2s 
								// (+/- 0.1s repeatability)

								*iters_ptr = p2;

								// This works suprisingly well- degradation is 
								// really only noticeable at high frequency 
								// transitions (e.g. with striped palettes).
								// Maybe average the mags at the 4 offsets to 
								// make it better

								// The mag store causes about a 7.5% slowdown 
								// on zoomtest.
								MAG(m, iters_ptr) = MAG(m, &iters_ptr[offs2]);

								// this adds no measureable overhead
								points_guessed++;
							}
							else
							{
								// Load RE coordinate from the array
								ps_ptr->ab_in[0] = m->img_re[x];
								
								m->queue_point(m, ps_ptr, iters_ptr);
							}
							iters_ptr += inc;
							x += inc;
						}
					} while ((y += inc) <= yend);
				}
				// really should flush at the end of each wave, but any errors 
				// should have no visual effect
			} // end of wave loop
		}
		s++; // go to next stripe
	} // end of stripe loop

	// accumulate iters, for thread load balance measurement
	t->total_iters += ps_ptr->iterctr;

	t->points_guessed = points_guessed;

	// Up to 4 points could be left in the queue (or 8 for SSE). Queue 
	// non-diverging dummy points to flush them out. This is tricky. Be careful 
	// changing it... can cause corrupted pixel bugs. Turns out that 4 more 
	// points (8 for SSE) must always be queued. They could be stored to the 
	// dummy value if all points left in the queue still have max_iters 
	// remaining.

	ps_ptr->ab_in[0] = 0.0;
	ps_ptr->ab_in[1] = 0.0;

	// Add some extra logic here to get the exact iteration count (i.e, exclude 
	// dummy iterations). It's actually pretty tough to calculate
	for (i = m->precision == PRECISION_SINGLE ? 8 : 4; i--;)
		m->queue_point(m, ps_ptr, m->iter_data + m->image_size);

	// Thread 0 always runs in the master thread, so doesn't need to signal. 
	// Save overhead.
	if (t->thread_num)
		SetEvent(t->done_event);	// For other threads, signal master thread 
									// that we're done

	return 0;
}

/**
 * Check for precision loss- occurs if the two doubles (or converted floats)
 * in ptest are equal to each other. Returns PLOSS_FLOAT for float loss,
 * PLOSS_DOUBLE for double loss, etc. or 0 for no loss.
 *
 * New more conservative version demands that bits beyond the lsb should also
 * differ. If only the lsb differs, bound to get degradation during iteration.
 */

#define PLOSS_DOUBLE    2
#define PLOSS_FLOAT     1

int check_precision_loss(double* ptest)
{
	float f[2];
	int i0[2], i1[2];

	// Check double loss
	i0[0] = *((int*)&ptest[0]) & ~1;	// get low dword of 1st double; mask 
										// off lsb
	i0[1] = *((int*)&ptest[0] + 1);		// get high dword

	i1[0] = *((int*)&ptest[1]) & ~1;	// get low dword of 2nd double; mask 
										// off lsb
	i1[1] = *((int*)&ptest[1] + 1);		// get high dword

	if ((i0[1] == i1[1]) && (i0[0] == i1[0]))
		return PLOSS_DOUBLE | PLOSS_FLOAT;	// double loss is also float loss

	 // Check float loss
	f[0] = (float)ptest[0];				// convert doubles to floats
	f[1] = (float)ptest[1];

	i0[0] = *((int*)&f[0]) & ~1;		// get 1st float; mask off lsb
	i1[0] = *((int*)&f[1]) & ~1;		// get 2nd float; mask off lsb

	if (i0[0] == i1[0])
		return PLOSS_FLOAT;

	return 0;
}

/**
 * Calculate the real and imaginary arrays for the current rectangle, set 
 * precision/algorithm, and do other misc setup operations. Call before 
 * starting mandelbrot calculation.
 *
 * Panning is now tracked as an offset from re/im, so the current rectangle is 
 * offset from re/im by (pan_xoffs + xstart / 2) and (pan_yoffs + ystart / 2).
 *
 * sms
 */
void man_setup(man_calc_struct* m, int xstart, int xend, int ystart, int yend)
{
	int i, x, y, xsize, ysize, ploss;
	long long step;
	unsigned int queue_init, flags;
	man_pointstruct* ps_ptr;

	m->max_iters &= ~1;	// make max iters even- required by optimized alg
	xsize = m->xsize;
	ysize = m->ysize;
	flags = m->flags;

	// Make re/im arrays, to avoid doing xsize * ysize flops in the main loop.
	// Also check for precision loss (two consecutive values equal or differing 
	// only in the lsb)

	// Cut overhead by only going from start to end - not whole image size
	// Need to do more than 1 point to detect precision loss... otherwise auto 
	// mode won't work.
	// To be safe use at least 4 (allocate 4 extra).

	// Updated to use offsets (pan_xoffs and pan_yoffs) for panning, rather 
	// than updating re/im on every pan. These are added in below. See comments 
	// at top (bug fix)

	// If saving, don't check precision loss, and don't recalculate the re 
	// array after the first row.

	if (!(flags & FLAG_IS_SAVE))
	{
		xend += 4;  // only need these for non-save (precision loss checking)
		yend += 4;
	}

	// this flag should be 1 for main calculation
	if (flags & FLAG_CALC_RE_ARRAY)
	{
		ploss = 0;
		x = xstart;
		step = -(xsize >> 1) + xstart + m->pan_xoffs;
		do
		{
			m->img_re[x] = m->re + get_re_im_offs(m, step++);
			if (!(flags & FLAG_IS_SAVE) && x > xstart)
				ploss |= check_precision_loss(&m->img_re[x - 1]);
		} while (++x <= xend);
	}

	step = -(ysize >> 1) + ystart + m->pan_yoffs;
	y = ystart;
	do
	{
		m->img_im[y] = m->im - get_re_im_offs(m, step++);
		if (!(flags & FLAG_IS_SAVE) && y > ystart)
			ploss |= check_precision_loss(&m->img_im[y - 1]);
	} while (++y <= yend);

	if (!(flags & FLAG_IS_SAVE)) // only do auto precision if not saving
	{
		precision_loss = 0;

		// Set precision loss flag. If in auto precision mode, set single or 
		// double calculation precision based on loss detection.
		switch (m->precision)
		{
		case PRECISION_AUTO:
			m->precision = PRECISION_SINGLE;
			if (ploss & PLOSS_FLOAT)
				m->precision = PRECISION_DOUBLE; // deliberate fallthrough
		case PRECISION_DOUBLE:
			if (ploss & PLOSS_DOUBLE)
				precision_loss = 1;
			break;
		case PRECISION_SINGLE:
			if (ploss & PLOSS_FLOAT)
				precision_loss = 1;
			break;
		default:
			// should never get here (x87 is suppressed until implemented)
			break;
		}
	}

	// Set iteration and queue_point function pointers and initialize queues

	// Alg will always be C if no sse support.
	// Should change algorithm in dialog box if it's reset to C here.
	if ((m->alg & ALG_C) || 
		(sse_support < 2 && m->precision == PRECISION_DOUBLE))
	{
		m->queue_point = queue_point_c;
		m->mandel_iterate = iterate_c; // Unoptimized C algorithm
	}
	else
	{
		if (m->precision == PRECISION_DOUBLE)
		{
			queue_init = (QUEUE_FULL << 12) | 
				(3 << 9) | (2 << 6) | (1 << 3) | 0;

			m->queue_point = queue_4point_sse2;
			m->mandel_iterate = (m->alg & ALG_INTEL) 
				? iterate_intel_sse2 : iterate_amd_sse2;
		}
		else
		{
			queue_init = (QUEUE_FULL << 24) | 
				(7 << 21) | (6 << 18) | (5 << 15) |
				(4 << 12) | (3 << 9) | (2 << 6) | (1 << 3) | 0;

			m->queue_point = queue_8point_sse;
			m->mandel_iterate = (m->alg & ALG_INTEL) 
				? iterate_intel_sse : iterate_amd_sse;
		}
	}

	// Set pointstruct initial values
	for (i = 0; i < num_threads; i++)
	{
		ps_ptr = m->thread_states[i].ps_ptr;
		ps_ptr->queue_status = queue_init;
		ps_ptr->cur_max_iters = m->max_iters;
		ps_ptr->iterctr = 0;
	}
}

/**
 * man_calculate() splits the calculation up into multiple threads, each 
 * calling the man_calculate_threaded() function.
 *
 * Current alg: divide the calculation rectangle into N stripes per thread (N
 * depends on the number of threads). More stripes help load balancing by 
 * making it unlikely that the stripes will have wildly different iteration 
 * counts. But too many stripes cause a slowdown due to excess overhead.
 *
 * Example, image rectangle and stripes for two threads:
 *
 *       1 Stripe            2 Stripes         etc...
 *
 *  +----------------+   +----------------+
 *  |                |   |   Thread 0     |
 *  |   Thread 0     |   |----------------|
 *  |                |   |   Thread 1     |
 *  +----------------|   |----------------|
 *  |                |   |   Thread 0     |
 *  |   Thread 1     |   |----------------|
 *  |                |   |   Thread 1     |
 *  +----------------+   +----------------+
 *
 *  The rectangle is sometimes divided horizontally:
 *
 *  +-------------------------------+
 *  |  T0   |  T1   |  T0   |  T1   |
 *  +-------------------------------+
 *
 * This is necessary so that 1-pixel high rectangles (as often found in 
 * panning) can still be divided. Horizontal stripes (vertical division) would 
 * be preferred because x is done in the inner loop. Also memory access is 
 * better on horizontal stripes. Arbitrarily decide to do horizontal division 
 * only if the stripe height is < some constant (say 8, = 2x fast alg cell 
 * size).
 *
 * Returns the time taken to do the calculation.
 *
 * smc
 */
double 
man_calculate(man_calc_struct* m, int xstart, int xend, int ystart, int yend)
{
	TIME_UNIT start_time;
	double iteration_time;
	int i, xsize, ysize, step, thread_ind, stripe_ind, num_stripes, frac;
	int frac_step, this_step;
	stripe* s = NULL;

	all_recalculated = 0;
	
	// if need recalculation, recalculate all. No effect for saving
	if (status & STAT_NEED_RECALC)
	{
		xstart = 0;	// reset rectangle to full screen
		xend = m->xsize - 1;
		ystart = 0;
		yend = m->ysize - 1;
		status &= ~STAT_NEED_RECALC;
		all_recalculated = 1;
	}

	man_setup(m, xstart, xend, ystart, yend);

	xsize = xend - xstart + 1;
	ysize = yend - ystart + 1;

	// Get the number of stripes per thread based on number of threads (extract 
	// bits 3-0 for 1 thread, 7-4 for two threads, etc).
	num_stripes = (cfg_default.stripes_per_thread.c_val >> (num_threads_ind << 2)) & 0xF;

	// Need to check min/max here (couldn't be checked automatically by 
	// log-reading function)
	if (num_stripes < 1)
		num_stripes = 1;
	if (num_stripes > MAX_STRIPES)
		num_stripes = MAX_STRIPES;

	// Now multiply by num_threads to get total number of stripes for the image.
	num_stripes <<= num_threads_ind;

	// With pathologically small images, some threads may not calculate anything.
	for (i = 0; i < num_threads; i++)
		m->thread_states[i].num_stripes = 0;

	// Start at the last thread, so that thread 0 gets any leftovers at the end. 
	// Thread 0 is the master and doesn't suffer the overhead of being spawned, 
	// so it should get the extra work.

	thread_ind = num_threads - 1;
	stripe_ind = 0;

	// Divide along the y axis if stripe height is >= 8 (see above), or ysize 
	// is >= xsize

	if ((ysize >= (num_stripes << 3)) || (ysize >= xsize))
	{
		// step size (height of each stripe)
		if (!(step = ysize / num_stripes))
		{
			num_stripes = ysize;	// if more stripes than pixels of height,
			step = 1;				// limit stripes and threads
		}

		// Use fractional steps to get the threads as evenly balanced as 
		// possible. For each thread, the stripe height could either be step or 
		// step + 1 (they all get the same step as a group).
		//
		// Dual Opteron 280, Double, Fast, 4 threads, 4 stripes/thread, tune.log:
		// With fractional steps	:	4.788s, efficiency 97.8%
		// Without fractional steps	:	4.980s, efficiency 94.1%

		frac = frac_step = ysize - (num_stripes * step);
		this_step = step;

		for (i = 0; i < num_stripes; i++)
		{
			m->thread_states[thread_ind].num_stripes++;
			s = &m->thread_states[thread_ind].stripes[stripe_ind];
			s->xstart = xstart;
			s->xend = xend;
			s->ystart = ystart;
			s->yend = ystart + this_step - 1;
			ystart += this_step;

			// Next stripe goes to next thread. If it wraps, reset and 
			// increment each thread's stripe index.
			if (--thread_ind < 0)
			{
				thread_ind = num_threads - 1;
				stripe_ind++;

				// Now that each thread has a stripe, update the fraction, and 
				// if it wraps increase the stripe height for the next group 
				// by 1.
				this_step = step;
				if ((frac += frac_step) >= num_stripes)
				{
					frac -= num_stripes;
					this_step++;	// comment this out to compare efficiencies 
									// without fractional steps
				}
			}
		}
		// Give thread 0 any leftovers.
		s->yend = yend;
	}
	else // Similar code to above for dividing along the x axis
	{
		if (!(step = xsize / num_stripes))
		{
			num_stripes = xsize;
			step = 1;
		}
		frac = frac_step = xsize - (num_stripes * step);
		this_step = step;

		for (i = 0; i < num_stripes; i++)
		{
			m->thread_states[thread_ind].num_stripes++;
			s = &m->thread_states[thread_ind].stripes[stripe_ind];
			s->xstart = xstart;
			s->xend = xstart + this_step - 1;
			s->ystart = ystart;
			s->yend = yend;
			xstart += this_step;

			if (--thread_ind < 0)
			{
				thread_ind = num_threads - 1;
				stripe_ind++;
				this_step = step;
				if ((frac += frac_step) >= num_stripes)
				{
					frac -= num_stripes;
					this_step++;
				}
			}
		}
		s->xend = xend;
	}

	/* Run the threads on the stripes calculated above.
	 * These threading functions are slow.
	 * Benchmarks on an Athlon 64 4000+ 2.4 GHz:
	 * 
	 * functions (tested with 4 threads):
	 * T1: 4 CreateThread + WaitForMultipleObjects + 4 CloseHandle
	 * T2: 4 _beginthreadex + WaitForMultipleObjects + 4 CloseHandle
	 * T3: 4 QueueUserWorkItem + 4 SetEvent + WaitForMultipleObjects
	 * T4: 4 SetEvent
	 * T5: 4 Null (loop overhead + mandel function call only)
	 * 
	 *														Equivalent SSE2
	 * Functions					uS Clock	Cycles		iters (per core)
	 * ------------------------------------------------------------------------
	 * T1							173.00		415200		83040
	 * T2							224.00		537600		107520
	 * T3							7.60		18240		3648
	 * T4							1.00		2400		480
	 * T5							0.01		24			4
	 *
	 * The QueueUserWorkItem method is probably fast enough. The first two are 
	 * usually on par with the iteration time while panning, negating any 
	 * multi-core advantage (for panning).

	 * Don't call C library routines in threads. 4K stack size is more than 
	 * enough
	 */
	start_time = get_timer();

	// Using WT_EXECUTEINPERSISTENTTHREAD is 50% slower than without.
	// Leaving out WT_EXECUTELONGFUNCTION makes the initial run twice as slow, 
	// then same speed.

	// Use the master thread (here) to do some of the work. Queue any other 
	// threads. Saves some overhead, and doesn't spawn any new threads at all 
	// if there's only one thread.

	for (i = 1; i < num_threads; i++)
		QueueUserWorkItem(man_calculate_threaded, &m->thread_states[i],
			WT_EXECUTELONGFUNCTION | (MAX_QUEUE_THREADS << 16));

	// Could also queue thread 0 too, then have this master thread display the 
	// progress if the calculation is really slow...
	man_calculate_threaded(&m->thread_states[0]);

	if (num_threads > 1)
	{
#ifdef _WIN32
		WaitForMultipleObjects(num_threads - 1, 
			&m->thread_done_events[1], 
			TRUE, INFINITE); // wait till all threads are done
#endif
	}

	if (!(m->flags & FLAG_IS_SAVE)) // don't update these if doing save
	{
		iteration_time = get_seconds_elapsed(start_time);
		file_tot_time += iteration_time;
		return iteration_time;
	}
	else
		return 0.0;
}

/* ---------------------- Quadrant/panning functions ----------------------- */

/**
 * Swap the memory pointers and handles of two quadrants (e.g., upper left and 
 * upper right).
 * Doesn't modify other fields as these will all be changed afterwards.
 */
void swap_quadrants(quadrant* q1, quadrant* q2)
{
	unsigned int* tmp_data;
	HBITMAP tmp_handle;

	tmp_data = q1->bitmap_data;
	tmp_handle = q1->handle;

	q1->bitmap_data = q2->bitmap_data;
	q1->handle = q2->handle;

	q2->bitmap_data = tmp_data;
	q2->handle = tmp_handle;
}

/**
 * Reset the quadrants to the initial state: put the screen in the UL quadrant, 
 * set the blit size to the screen size, and set one update rectangle to the UL 
 * quadrant. Set all other quadrants inactive. Calling this will cause 
 * recalculation of the whole image. Call if this is desired, or if the screen 
 * size changes.
 */
void reset_quadrants(void)
{
	int xsize, ysize;
	quadrant* q;
	man_calc_struct* m = &main_man_calc_struct;

	xsize = m->xsize;
	ysize = m->ysize;

	q = &quad[UL];				// upper left
	q->status = QSTAT_DO_BLIT;	// needs blit
	q->src_xoffs = 0;			// src and dest offsets all 0
	q->src_yoffs = 0;
	q->dest_xoffs = 0;
	q->dest_yoffs = 0;
	q->blit_xsize = xsize;		// blit sizes = screen size
	q->blit_ysize = ysize;
	q->quad_rect.x[0] = 0;		// rectangle coordinates
	q->quad_rect.y[0] = 0;
	q->quad_rect.x[1] = xsize - 1;
	q->quad_rect.y[1] = ysize - 1;

	q = &quad[UR];				// upper right
	q->status = 0;				// inactive
	q->quad_rect.x[0] = xsize;	// rectangle coordinates
	q->quad_rect.y[0] = 0;
	q->quad_rect.x[1] = (xsize << 1) - 1;
	q->quad_rect.y[1] = ysize - 1;

	q = &quad[LL];				// lower left
	q->status = 0;				// inactive
	q->quad_rect.x[0] = 0;		// rectangle coordinates
	q->quad_rect.y[0] = ysize;
	q->quad_rect.x[1] = xsize - 1;
	q->quad_rect.y[1] = (ysize << 1) - 1;

	q = &quad[LR];				// lower right
	q->status = 0;				// inactive
	q->quad_rect.x[0] = xsize;	// rectangle coordinates
	q->quad_rect.y[0] = ysize;
	q->quad_rect.x[1] = (xsize << 1) - 1;
	q->quad_rect.y[1] = (ysize << 1) - 1;

	// Update rectangles
	update_rect[0].valid = 1;	// 1st update rect valid; equals whole quadrant
	update_rect[0].x[0] = 0;
	update_rect[0].x[1] = xsize - 1;
	update_rect[0].y[0] = 0;
	update_rect[0].y[1] = ysize - 1;
	update_rect[1].valid = 0;	// 2nd update rect invalid

	screen_xpos = 0;	// reset screen window to cover UL quadrant only
	screen_ypos = 0;
}

/**
 * Calculate the intersection of two rectangles. If they intersect, sets rdest 
 * to the intersection and returns 1, else returns 0.
 * For each rectangle, the coord at index 0 is <= the coord at index 1.
 */
int intersect_rect(rectangle* rdest, rectangle* r1, rectangle* r2)
{
	// Check if one lies outside the bounds of the other; return 0 if so
	if (r1->x[0] > r2->x[1] || r1->x[1] < r2->x[0] ||
		r1->y[0] > r2->y[1] || r1->y[1] < r2->y[0])
		return 0;

	// Get max/min coordinates to get intersection
	rdest->x[0] = r1->x[0] > r2->x[0] ? r1->x[0] : r2->x[0];
	rdest->x[1] = r1->x[1] < r2->x[1] ? r1->x[1] : r2->x[1];
	rdest->y[0] = r1->y[0] > r2->y[0] ? r1->y[0] : r2->y[0];
	rdest->y[1] = r1->y[1] < r2->y[1] ? r1->y[1] : r2->y[1];
	return 1;
}

/**
 * Iterate on the update rectangles, and palette-map the iteration data to the 
 * quadrants. Only used for main calculation, not while saving.
 *
 * smq
 */
void man_calculate_quadrants(void)
{
	rectangle r;
	int i, j, x, y;
	unsigned int* bmp_ptr, * iters_ptr;
	man_calc_struct* m = &main_man_calc_struct;

	iter_time = 0.0;

	// First calculate the update rectangles (up to 2).
	for (i = 0; i < 2; i++)
		if (update_rect[i].valid)
		{
			// To get position in (screen-mapped) image, subtract screen upper 
			// left coordinates, Rectangles will be at one of the screen edges 
			// (left, right, top, or bottom).
			// Could simplify this: determined solely by pan offs_x and offs_y

			// Iterate on the update rectangles
			iter_time += man_calculate
			(m, update_rect[i].x[0] - screen_xpos,  // xstart
				update_rect[i].x[1] - screen_xpos,  // xend
				update_rect[i].y[0] - screen_ypos,  // ystart
				update_rect[i].y[1] - screen_ypos); // yend
		}

	// Now palette-map the update rectangles into their quadrants. Each 
	// rectangle can occupy 1-4 quadrants.

	for (i = 0; i < 4; i++)
		for (j = 0; j < 2; j++)
			if (update_rect[j].valid)
				if (intersect_rect(&r, &quad[i].quad_rect, &update_rect[j]))
				{
					// Subtract upper left coordinates of this quadrant from 
					// upper left coords of intersection rect to get the x, y 
					// offset of the bitmap data in this quadrant.

					x = r.x[0] - quad[i].quad_rect.x[0];
					y = r.y[0] - quad[i].quad_rect.y[0];

					// get pointer to bitmap data
					bmp_ptr = quad[i].bitmap_data + y * m->xsize + x;

					// Subtract upper left coords of screen pos from upper left 
					// coords of intersection rect to get iter data offset

					x = r.x[0] - screen_xpos;
					y = r.y[0] - screen_ypos;

					// get pointer to iter data to be mapped
					iters_ptr = m->iter_data + y * m->iter_data_line_size + x;

					// Xsize, ysize = rectangle edge lengths
					apply_palette(m, bmp_ptr, iters_ptr, 
						r.x[1] - r.x[0] + 1, 
						r.y[1] - r.y[0] + 1, num_threads);
				}
}

/**
 * Pan the image by offs_x and offs_y. Sets iter_time (from iteration functions)
 * to the time to do the iteration only.
 */
void pan_image(int offs_x, int offs_y)
{
	int xsize, ysize, swap_x, swap_y, tmp;
	quadrant* q;
	rectangle* u;
	man_calc_struct* m = &main_man_calc_struct;

	if (offs_x | offs_y) // Recalculate only if the image moved
	{
		// Update pan offsets
		m->pan_xoffs -= offs_x;	// maybe invert offs_x, offs_y signs later
		m->pan_yoffs -= offs_y;

		xsize = m->xsize;
		ysize = m->ysize;

		// See algorithm explanation above
		screen_xpos -= offs_x;  // update screen pos
		screen_ypos -= offs_y;

		// Renormalize screen coordinates and swap quadrants if necessary
		swap_x = swap_y = 0;

		if (screen_xpos < 0)
		{
			screen_xpos += xsize;
			swap_x = 1;
		}
		if (screen_xpos > xsize)
		{
			screen_xpos -= xsize;
			swap_x = 1;
		}
		if (screen_ypos < 0)
		{
			screen_ypos += ysize;
			swap_y = 1;
		}
		if (screen_ypos > ysize)
		{
			screen_ypos -= ysize;
			swap_y = 1;
		}
		if (swap_x)
		{
			swap_quadrants(&quad[UL], &quad[UR]);
			swap_quadrants(&quad[LL], &quad[LR]);
		}
		if (swap_y)
		{
			swap_quadrants(&quad[UL], &quad[LL]);
			swap_quadrants(&quad[UR], &quad[LR]);
		}

		/* Get the update rectangles (there can be either 1 or 2): determined 
		   by difference between screen pos and previous screen pos (offs_x and 
		   offs_y). Use quadrant coordinates (0,0 to 2 * xsize -1, 2 * ysize -1). 
		   Rectangles will be split at quadrant boundaries during palette 
		   mapping. If either of the sizes is 0, the rectangle will be ignored 
		   (not calculated or palette mapped). */

		// Check for vertical rectangles (have an x offset)
		u = &update_rect[1];
		u->valid = offs_x;                  // valid if any nonzero offset
		u->y[0] = screen_ypos;              // vertical rectangles are
		u->y[1] = screen_ypos + ysize - 1;  // full height of screen

		if (offs_x > 0)
		{  // vertical rect at left
			u->x[0] = screen_xpos;
			u->x[1] = screen_xpos + offs_x - 1;
		}
		if (offs_x < 0)
		{  // vertical rect at right
			u->x[0] = screen_xpos + xsize + offs_x;
			u->x[1] = screen_xpos + xsize - 1;
		}

		// Check for horizontal rectangles (have a y offset).
		// Optimize out the corner intersection with any vertical rectangles so 
		// we don't calculate it twice. Clip it off the vertical rect. The 
		// intersection could be a large area for drag pans.

		u = &update_rect[0];
		u->valid = offs_y;                  // valid if any nonzero offset
		u->x[0] = screen_xpos;              // horizontal rectangles are
		u->x[1] = screen_xpos + xsize - 1;  // full width of screen

		if (offs_y > 0)
		{  // horizontal rect at top
			u->y[0] = screen_ypos;
			u->y[1] = tmp = screen_ypos + offs_y - 1;
			
			// clip off corner intersection from any vertical rect
			u[1].y[0] = tmp + 1;
		}
		if (offs_y < 0)
		{  // horizontal rect at bottom
			u->y[0] = tmp = screen_ypos + ysize + offs_y;
			u->y[1] = screen_ypos + ysize - 1;
			
			// clip off corner intersection from any vertical rect
			u[1].y[1] = tmp - 1;
		}

		// Get the blit rectangles from the screen position (screen_xpos, 
		// screen_ypos = screen upper left corner). Screen coordinates always 
		// range from 0 to xsize and 0 to ysize inclusive here.
		// Refer to diagram.

		quad[UL].status = 0; // Default: all inactive
		quad[UR].status = 0;
		quad[LL].status = 0;
		quad[LR].status = 0;

		// Check if UL has a blit rectangle
		if (screen_xpos < xsize && screen_ypos < ysize)
		{
			q = &quad[UL];
			q->status = QSTAT_DO_BLIT; // need a blit
			q->dest_xoffs = 0; // always blits to screen upper left corner
			q->dest_yoffs = 0;
			q->src_xoffs = screen_xpos;
			q->src_yoffs = screen_ypos;
			q->blit_xsize = xsize - screen_xpos;
			q->blit_ysize = ysize - screen_ypos;
		}
		// Check if UR has a blit rectangle
		if (screen_xpos > 0 && screen_ypos < ysize)
		{
			q = &quad[UR];
			q->status = QSTAT_DO_BLIT; // need a blit
			q->dest_xoffs = xsize - screen_xpos;
			q->dest_yoffs = 0;	// always blits to screen upper edge
			q->src_xoffs = 0;	// always blits from bitmap left edge
			q->src_yoffs = screen_ypos;
			q->blit_xsize = screen_xpos;
			q->blit_ysize = ysize - screen_ypos;
		}
		// Check if LL has a blit rectangle
		if (screen_xpos < xsize && screen_ypos > 0)
		{
			q = &quad[LL];
			q->status = QSTAT_DO_BLIT; // need a blit
			q->dest_xoffs = 0;         // always blits to screen left edge
			q->dest_yoffs = ysize - screen_ypos;
			q->src_xoffs = screen_xpos;
			q->src_yoffs = 0;          // always blits from bitmap top edge
			q->blit_xsize = xsize - screen_xpos;
			q->blit_ysize = screen_ypos;
		}
		// Check if LR has a blit rectangle
		if (screen_xpos > 0 && screen_ypos > 0)
		{
			q = &quad[LR];
			q->status = QSTAT_DO_BLIT; // need a blit
			q->dest_xoffs = xsize - screen_xpos;
			q->dest_yoffs = ysize - screen_ypos;
			q->src_xoffs = 0; // always blits from bitmap upper left corner
			q->src_yoffs = 0;
			q->blit_xsize = screen_xpos;
			q->blit_ysize = screen_ypos;
		}

		status |= STAT_RECALC_FOR_PALETTE;

		// Calculate and palette map the new update rectangles
		do_man_calculate(0);
	}
}

/**
 * Create the text in the INFO area of the main dialog
 * If update_iters_sec is 0, won't update iters/sec and gflops (use during
 * panning, when calculations will be inaccurate due to short calculation times)
 */
char* get_image_info(int update_iters_sec)
{
	static char s[1024 + 32 * MAX_THREADS];
	static char iters_str[256];
	static unsigned long long ictr = 0;
	static double guessed_pct = 0.0;
	static double miters_s = 0.0;	// mega iterations/sec
	static double avg_iters = 0.0;	// average iterations per pixel

	unsigned long long ictr_raw;
	unsigned long long ictr_total_raw;
	double cur_pct, max_cur_pct, tot_pct, max_tot_pct;
	int i, points_guessed;
	thread_state* t;
	char tmp[256];
	man_calc_struct* m = &main_man_calc_struct;

	ictr_raw = 0;
	ictr_total_raw = 0;
	points_guessed = 0;
	for (i = 0; i < num_threads; i++)
	{
		t = &m->thread_states[i];
		points_guessed += t->points_guessed;
		ictr_raw += t->ps_ptr->iterctr;	// N iterations per tick
		ictr_total_raw += t->total_iters;
	}

	if (update_iters_sec)
	{
		// For the C versions, each tick is 1 iteration.
		// For the SSE2 ASM versions, each tick is 4 iterations.
		// For the SSE ASM versions, each tick is 8 iterations.

		ictr = ictr_raw;
		if (!(m->alg & ALG_C))
		{
			if (m->precision == PRECISION_DOUBLE && sse_support >= 2)
				ictr <<= 2;
			if (m->precision == PRECISION_SINGLE && sse_support >= 1)
				ictr <<= 3;
		}

		// Prevent division by 0. If the time is in this neighborhood the 
		// iters/sec won't be accurate anyway.
		if (iter_time < 0.001)
			iter_time = 0.001;

		miters_s = (double)ictr * 1e-6 / iter_time;
		avg_iters = (double)ictr / (double)m->image_size;

		guessed_pct = 100.0 * (double)points_guessed / (double)m->image_size;

		// Since one flop is optimized out per 18 flops in the ASM versions,
		// factor should really be 8.5 for those. But actually does 9 
		// "effective" flops per iteration

		sprintf_s(iters_str, sizeof(iters_str), "%-4.4gM (%-.2f GFlops)", 
			miters_s, miters_s * 9.0 * 1e-3);
	}

	sprintf_s(s, sizeof(s),   // Microsoft wants secure version
		"Real\t%-16.16lf\r\n"
		"Imag\t%-16.16lf\r\n"
		"Mag\t%-16lf\r\n"
		"\r\n"
		"Size\t%u x %u\r\n"
		"Time\t%-3.3fs\r\n"
		"Iters/s\t%s\r\n"

		"\r\n" // Everything here and below is "hidden" - must scroll down to see

		"Avg iters/pixel\t%-.1lf\r\n"
		"Points guessed\t%-.1lf%%\r\n"
		"Total iters\t%-.0lf\r\n",

		// With new panning method, need to get actual screen centerpoint using 
		// pan offsets
		m->re + get_re_im_offs(m, m->pan_xoffs),
		m->im - get_re_im_offs(m, m->pan_yoffs),
		
		// Miters/s string created above
		m->mag, m->xsize, m->ysize, iter_time, iters_str, avg_iters, 
		guessed_pct, (double)ictr
	);

	// Get each thread's percentage of the total load, to check balance.

	sprintf_s(tmp, sizeof(tmp), "\r\nThread load %%\tCur    Total\r\n");
	strcat_s(s, sizeof(s), tmp);

	max_cur_pct = 0.0;
	max_tot_pct = 0.0;
	for (i = 0; i < num_threads; i++)
	{
		t = &m->thread_states[i];
		sprintf_s(tmp, sizeof(tmp),
			"Thread %d\t%#3.3g   %#3.3g\r\n", i,
			cur_pct = (double)t->ps_ptr->iterctr / (double)ictr_raw * 100.0,
			tot_pct = (double)t->total_iters / (double)ictr_total_raw * 100.0);
		strcat_s(s, sizeof(s), tmp);

		if (cur_pct > max_cur_pct)
			max_cur_pct = cur_pct;
		if (tot_pct > max_tot_pct)
			max_tot_pct = tot_pct;
	}

	// Figure of merit: percentage of best possible speed, which occurs with
	// perfect thread load balancing.
	sprintf_s(tmp, sizeof(tmp),
		"Efficiency %%\t%#3.3g   %#3.3g\r\n"
		"\r\nTotal calc time\t%-.3lfs\r\n" // resets whenever new file is opened
		"\r\n(C) 2006-2008 Paul Gentieu",
		100.0 * 100.0 / (num_threads * max_cur_pct),
		100.0 * 100.0 / (num_threads * max_tot_pct),
		file_tot_time);
	strcat_s(s, sizeof(s), tmp);

	return s;
}

/**
 * Print frames/sec info and current algorithm indicator. Current algorithm can
 * temporarily change from Fast to Exact while panning, to improve performance.
 * See man_calculate_threaded.
 */
void print_fps_status_line(double fps, double avg_fps, double eff)
{
	char s[256];
	man_calc_struct* m = &main_man_calc_struct;

	// if currently saving, keep saving status line
	if (status & STAT_DOING_SAVE)
		return;

	// Interval and average frames/sec. Removed "AVG" so large frame rates 
	// don't get cut off
	sprintf_s(s, sizeof(s), "%c Fps %3.0f/%-3.0f", 
		m->cur_alg & ALG_EXACT ? 'E' : 'F', fps, avg_fps);
	SetWindowText(hwnd_status, s);

	// Iteration percentage: mandelbrot calculation time / total time
	sprintf_s(s, sizeof(s), "Iter %2.0f%%", eff);
	SetWindowText(hwnd_status2, s);
}

/**
 * Get frames/sec during one interval, average frames per sec, and iteration 
 * time percentage (mandelbrot calculation time / pan or zoom time). Assumes 
 * iter_time is set to the mandelbrot calculation time. Op_time should be the 
 * total time for the pan or zoom operation. Call once per frame.
 */
void update_benchmarks(double op_time, int update_iters_sec)
{
#define UPDATE_INTERVAL_TIME  0.25 // update about 4 times per sec

	interval_frames++;
	total_frames++;
	interval_time += op_time;	// Update interval and total operation times
	total_time += op_time;
	
	// Update time spent only on mandelbrot calculation
	calc_total_time += iter_time;
	
	// Update time spent only on mandelbrot calculation this interval
	calc_interval_time += iter_time;

	// Update status line every interval
	if (interval_time >= UPDATE_INTERVAL_TIME)
	{
		print_fps_status_line((double)interval_frames / interval_time,
			(double)total_frames / total_time, 
			100.0 * (double)calc_interval_time / interval_time); // interval iteration %
			/*100.0 * (double) calc_total_time / total_time);     // average iteration %*/

		// Update image info
		SetWindowText(hwnd_info, get_image_info(update_iters_sec));
		
		interval_frames = 0;
		interval_time = 0.0;
		calc_interval_time = 0.0;
	}
}

// Some magic constants
#define PAN_STEP_DIV		150000.0
#define OVERHEAD_FACTOR		100000

/**
 * These adjust the pan filter constant based on image size (filter constant 
 * now comes from config file- adds acceleration and deceleration to movements).
 */
#define PFC_SLOPE_FACTOR   (1600.0 * 1140.0 - 700.0 * 700.0)
#define PFC_OFFS_FACTOR    (700.0 * 700.0)

static double cur_pan_xstep = 0.0, cur_pan_ystep = 0.0;
static double pan_xstep_accum = 0.0, pan_ystep_accum = 0.0;

#define PAN_KEY (KEY_RIGHT | KEY_LEFT | KEY_UP | KEY_DOWN) // any pan key

/* Reset the pan filter state and accumulators */
void reset_pan_state(void)
{
	cur_pan_xstep = 0.0;
	cur_pan_ystep = 0.0;
	pan_xstep_accum = 0.0;
	pan_ystep_accum = 0.0;
}

/* This calculates the x and y steps for panning based on keys pressed.
 * Returns 0 if it didn't do anything (due to loss of focus or whatever),
 * else 1.
 *
 * With very slow panning, the steps will be 0 most of the time, which will
 * cause the main code to automatically give up CPU time.
 *
 * Later make a more sophisticated algorithm that adjusts in real time,
 * based on the image characteristics and CPU speed.
 *
 * Returns a 1-clock pulse on a transition from active to stopped (keys haven't
 * been pressed for awhile, and all filters are basically cleared out).
 */
int get_pan_steps(int* xstep, int* ystep, int set_pan_key)
{
	int key, pkey, xstep_int, ystep_int, pulse;
	static int key_lock = 0, wait_release;
	double pan_step, xs, ys, pan_step_scale, pan_filter_const, tmp;
	double pfcmin, pfcmax, pfc_slope, pfc_offs;
	static int stopped_counter = 0;
	static int stopped = 0;
	man_calc_struct* m = &main_man_calc_struct;

#define STOPPED_COUNTER_MAX 25

	if (xstep == NULL || ystep == NULL) // null pointers mean set pan key
	{
		key_lock = set_pan_key;
		return 0;
	}

	// Without this, will pan even if other applications have the focus... 
	// still does if it recaptures focus due to the cursor moving over the main 
	// window.

	if (GetFocus() != hwnd_main)
	{
		// Stay active even if we don't have the focus, but don't accept any 
		// keys.
		// Necessary for clearing/updating the pan filters, and for keeping pan 
		// lock going. Maybe go into low-priority mode too?
		key = pkey = 0;
	}
	else
		key = pkey = get_keys_pressed();

	// No need to do all the code below if we're stopped and no pan keys were 
	// pressed.
	if (!(key & (PAN_KEY | KEY_CTRL)) && !key_lock)
	{
		if (stopped)
			return 0;
	}
	else // a pan key was pressed
		stopped = 0;

	// Do pan lock with the shift key. If pan is locked, keep panning in the
	// current direction. CTRL can toggle speed. Fun...
	if (key_lock)
	{
		if (pkey & (PAN_KEY | KEY_CTRL))
		{
			if (!wait_release)
				if (pkey & KEY_CTRL)
				{
					key_lock ^= KEY_CTRL;
					wait_release = 1;
				}
				else
					// preserve CTRL (fast/slow) state
					key_lock = (key_lock & ~PAN_KEY) | (pkey & PAN_KEY);
		}
		else
		{
			// leave lock mode if shift pressed without a pan key
			if (pkey & KEY_SHIFT)
				key_lock = 0;
			wait_release = 0;
		}
		key = key_lock;
	}
	else // not in pan lock mode: see if we need to go into it
		if ((pkey & KEY_SHIFT) && (pkey & PAN_KEY))
		{
			key_lock = pkey & (PAN_KEY | KEY_CTRL);
			wait_release = 1;
		}

	// Get the step scale from the pan rate slider
	pan_step_scale = pan_step_scales[cfg_default.pan_rate.c_val];

	// Adjust the step based on the image size, to try to compensate for 
	// size-based variable frame rates
	pan_step = pan_step_scale * (double)
		(m->image_size + OVERHEAD_FACTOR) * (1.0 / PAN_STEP_DIV);

	// Pan filter settings now come from config file. The 0.0001 scales large 
	// config file integers to the values needed here.
	pfcmin = 0.0001 * (double)cfg_default.pfcmin.c_val;
	pfcmax = 0.0001 * (double)cfg_default.pfcmax.c_val;
	pfc_slope = (pfcmax - pfcmin) * (1.0 / PFC_SLOPE_FACTOR);
	pfc_offs = pfcmin - pfc_slope * PFC_OFFS_FACTOR;

	pan_filter_const = (double)m->image_size * pfc_slope + pfc_offs;
	if (pan_filter_const < pfcmin)
		pan_filter_const = pfcmin;

	// Go to fast mode while CTRL key is held down
	if (key & KEY_CTRL)
		pan_step *= 4.0;	// arbitrary speedup factor: 4 seems to work pretty 
							// well

	 // Get x and y steps
	xs = key & KEY_RIGHT ? -pan_step : key & KEY_LEFT ? pan_step : 0.0;
	ys = key & KEY_DOWN ? -pan_step : key & KEY_UP ? pan_step : 0.0;

	// Filter the pan movements
	cur_pan_xstep = xs * pan_filter_const + (tmp = 1.0 - pan_filter_const) * cur_pan_xstep;
	cur_pan_ystep = ys * pan_filter_const + tmp * cur_pan_ystep;

	// Accumulate fractional steps
	pan_xstep_accum += cur_pan_xstep;
	pan_ystep_accum += cur_pan_ystep;

	// Round up/down based on sign; convert to int
	xstep_int = (int)(pan_xstep_accum + ((pan_xstep_accum < 0.0) ? -0.5 : 0.5));
	ystep_int = (int)(pan_ystep_accum + ((pan_ystep_accum < 0.0) ? -0.5 : 0.5));

	// Subtract integer part
	pan_xstep_accum -= (double)xstep_int;
	pan_ystep_accum -= (double)ystep_int;

	// Set integer steps
	*xstep = xstep_int;
	*ystep = ystep_int;

	// Detect when a pan stopped
	pulse = 0;
	// Reset the stopped counter on any activity
	if (key | xstep_int | ystep_int)
		stopped_counter = STOPPED_COUNTER_MAX;
	else
		// else decrement and pulse on transition to 0
		if (stopped_counter && !--stopped_counter)
		{
			stopped = 1;
			
			// Clear out old data. Without this, there can be little artifacts 
			// when restarting.
			reset_pan_state();
			pulse = 1;
		}

	return pulse;
}

/**
 * Pan the image using the keyboard. Super cool...
 *
 * Returns 1 if it did a pan, else 0 (if idle).
 */
int do_panning(void)
{
	int xstep = 0, ystep = 0;
	static TIME_UNIT start_time;
	static double pan_time = -1.0;

	// Update coords, xstep, and ystep based on keys pressed. Returns a 1-clock 
	// pulse when panning transitioned from active to stopped.

	if (get_pan_steps(&xstep, &ystep, 0))
	{
		SetWindowText(hwnd_info, get_image_info(0)); // update image info
		pan_time = -1.0;	// restart timing next time
		return 0;
	}

	// The Sleep value in the main loop can affect the timing when panning has 
	// almost decelerated to a stop (when the x and y steps are both zero, but 
	// may be 1 next cycle). In this case the Sleep time substitutes for the 
	// iteration + bitblit time. Want these to be as close as possible for 
	// smooth stops. Using a sleep value of 1 caused some discontinuities 
	// (sometimes doesn't sleep at all?) whereas 2 seems pretty good.
	// Later, probably want to do dummy iteration instead.

	if (xstep | ystep)
	{
		if (pan_time < 0.0) // if < 0, we need to restart timer
			start_time = get_timer();

		pan_image(xstep, ystep);

		// get time since last screen update
		pan_time = get_seconds_elapsed(start_time);

		// skip update if the whole image was recalculated (due to size change, 
		// etc). messes up the averages.
		if (!all_recalculated)
			update_benchmarks(pan_time, 0);	// update the fps, average fps, and 
											// iteration %
		start_time = get_timer();	// restart timer for next update

		return 1;
	}
	return 0;
}

// Get re/im coordinates at the mouse position (mx, my), for realtime zoom
void get_mouse_re_im(int mx, int my)
{
	man_calc_struct* m = &main_man_calc_struct;

	mx -= (m->xsize >> 1); // Get offset from image center
	my -= (m->ysize >> 1);

	mouse_re = m->re + get_re_im_offs(m, mx);
	mouse_im = m->im - get_re_im_offs(m, my);
}

/**
 * Do realtime zooming. Has two modes:
 *
 * Mode #1: mouse-controlled zooming. left button = zoom in, right button = 
 * zoom out.
 * Keeps the point under the mouse at a constant position on the screen.
 *
 * Mode #2: do a realtime zoom to the current image. This is done when the zoom
 * button is clicked.
 */
int do_zooming(void)
{
	TIME_UNIT start_time;
	int mx, my, done = 0;
	double step;
	man_calc_struct* m = &main_man_calc_struct;

	// If a panning key is pressed, temporarily exit to do the pan, then resume 
	// any zooming (but abort zooming started with the zoom button). Need to 
	// improve this...
	// Problem is that zoom frame rate is vastly lower than pan rate.

	if (get_keys_pressed() & PAN_KEY)
	{
		if (do_rtzoom)
			if (!(do_rtzoom & RTZOOM_WITH_BUTTON))
				prev_do_rtzoom = do_rtzoom;
			else
				prev_do_rtzoom = 0;
		do_rtzoom = 0;
	}
	else
	{
		if (prev_do_rtzoom)
		{
			// update re/im from any pan offsets and reset offsets
			update_re_im(m, m->pan_xoffs, m->pan_yoffs);

			// update mouse coords after any pan
			get_mouse_re_im(mouse_x[1], mouse_y[1]);
			reset_pan_state();	// reset pan so we don't get extra movement
			do_rtzoom = prev_do_rtzoom;	// after stopping zoom
		}
	}
	if (!do_rtzoom) // Return if not zooming
		return 0;

	// update re/im from any pan offsets and reset offsets
	update_re_im(m, m->pan_xoffs, m->pan_yoffs);

	step = rtzoom_mag_steps[cfg_default.zoom_rate.c_val];

	start_time = get_timer();

	if (do_rtzoom & RTZOOM_IN)
		m->mag *= step;
	else
	{
		m->mag /= step;
		if (m->mag < MAG_MIN)
			m->mag = MAG_MIN;
	}
	if (!(do_rtzoom & RTZOOM_WITH_BUTTON)) // if zooming using the mouse
	{
		// Set the new image center re/im to keep the position at mouse[1]
		// at the same point on the screen.

		mx = mouse_x[1] - (m->xsize >> 1); // Get offset from image center
		my = mouse_y[1] - (m->ysize >> 1);

		m->re = mouse_re - get_re_im_offs(m, mx);
		m->im = mouse_im + get_re_im_offs(m, my);
	}
	else // if zooming using the button, stop when we hit the start mag
		if (m->mag > zoom_start_mag)
		{
			m->mag = zoom_start_mag;
			done = 1;	// setting do_rtzoom 0 here wipes out fps numbers after 
						// button zoom is done
		}

	do_man_calculate(1);

	update_benchmarks(get_seconds_elapsed(start_time), 1);

	if (done) // If we just finished button zoom, update info
	{
		do_rtzoom = 0;
		// use for benchmarking
		file_tot_time = get_seconds_elapsed(zoom_start_time);

		SetWindowText(hwnd_info, get_image_info(1)); // Update all info
	}
	return 1;
}

/* Simple function for recalculating after a window size change (if enabled) */
int do_recalc(void)
{
	if (status & STAT_RECALC_IMMEDIATELY)
	{
		do_man_calculate(1);
		status &= ~STAT_RECALC_IMMEDIATELY;
	}
	return 0;
}

/**
 * Initialize values that never change. Call once at the beginning of the 
 * program. 
 */
void man_init(void)
{
	int i, j;
	HANDLE e;
	man_calc_struct* m;
	man_pointstruct* ps_ptr;

	// Initialize both main and save calculation structures
	for (i = 0; i < 2; i++)
	{
		m = i ? &save_man_calc_struct : &main_man_calc_struct;

		m->flags         = i 
			? FLAG_IS_SAVE | FLAG_CALC_RE_ARRAY : FLAG_CALC_RE_ARRAY;
		m->palette       = DEFAULT_PAL;
		m->rendering_alg = cfg_default.options.c_val & OPT_NORMALIZED 
			? RALG_NORMALIZED : RALG_STANDARD;
		m->precision     = PRECISION_AUTO;
		m->mag           = HOME_MAG;
		m->max_iters     = HOME_MAX_ITERS;

		// Initialize the thread state structures
		for (j = 0; j < MAX_THREADS; j++)
		{
			m->thread_states[j].thread_num = j;
			m->thread_states[j].calc_struct = m;

			// Create an auto-reset done event for each thread. 
			// The thread sets it when done with a calculation.
			e = CreateEvent(NULL, FALSE, FALSE, NULL);
			m->thread_states[j].done_event = e;
			m->thread_done_events[j] = e;

			// Init each thread's point structure
			m->thread_states[j].ps_ptr = ps_ptr = &m->pointstruct_array[j];

			// Init 64-bit double and 32-bit float fields with divergence 
			// radius and constant 2.0
			ps_ptr->two_d[1] = ps_ptr->two_d[0] = 2.0;
			ps_ptr->two_f[3] = ps_ptr->two_f[2] = ps_ptr->two_f[1] = 
				ps_ptr->two_f[0] = 2.0;

			ps_ptr->rad_d[1] = ps_ptr->rad_d[0] = DIVERGED_THRESH;
			ps_ptr->rad_f[3] = ps_ptr->rad_f[2] = ps_ptr->rad_f[1] = 
				ps_ptr->rad_f[0] = DIVERGED_THRESH;
		}
	}
}

/* ------------------------ GUI / misc functions --------------------------- */

/**
 * Detect whether the CPU supports SSE2 and conditional move instructions; used 
 * to set algorithms. Also detect the number of cores.
 */

#define FEATURE_SSE		0x02000000
#define FEATURE_SSE2	0x04000000
#define FEATURE_CMOV	0x00008000

void get_cpu_info(void)
{
	unsigned int vendor[4];
	unsigned int features;
	SYSTEM_INFO info;
	man_calc_struct* m = &main_man_calc_struct;

	vendor[3] = 0;

	__asm
	{
		// If the (ancient) CPU doesn't support the CPUID instruction, we would 
		// never get here...

		xor eax, eax	// get vendor
		cpuid
		mov vendor, ebx
		mov vendor + 4, edx
		mov vendor + 8, ecx
		mov eax, 1		// get features
		cpuid
		mov features, edx
	}

	// Use vendor to select default algorithm
	if (!strcmp((const char*)vendor, "AuthenticAMD"))
		m->alg = ALG_FAST_ASM_AMD;
	else
		m->alg = ALG_FAST_ASM_INTEL;

	// Set exact alg if configured by quickman.cfg options field
	if (cfg_default.options.c_val & OPT_EXACT_ALG)
		m->alg |= ALG_EXACT;

	sse_support = 0;
	if ((features & (FEATURE_SSE | FEATURE_CMOV)) == (FEATURE_SSE | FEATURE_CMOV))
		sse_support = 1;
	if ((features & (FEATURE_SSE2 | FEATURE_CMOV)) == (FEATURE_SSE2 | FEATURE_CMOV))
		sse_support = 2;

	if (sse_support < 2)
	{
		MessageBox(NULL, 
			"Your (obsolete) CPU does not support SSE2 instructions.\r\n"
			"Performance will be suboptimal.", "Warning", 
			MB_OK | MB_ICONSTOP | MB_TASKMODAL);

		// Ok to stay in auto precision mode with only sse support- will switch
		// to C algorithm for double. If no sse support, no choice but to use C
		if (!sse_support)
			m->alg = ALG_FAST_C;
	}

	// Set the default number of threads to the number of cores. Does this 
	// count a hyperthreading single core as more than one core? Should ignore 
	// these as hyperthreading won't help.
	GetSystemInfo(&info);
	num_threads = info.dwNumberOfProcessors;

	// Convert number of threads (cores) to a selection index for the dropdown 
	// box

	// Get log2(num_threads)
	for (num_threads_ind = 0; num_threads_ind <= MAX_THREADS_IND; num_threads_ind++)
		if ((1 << num_threads_ind) >= num_threads)
			break;

	num_threads = 1 << num_threads_ind;
}

/**
 * Allocate all the memory needed by the calculation engine. This needs to be 
 * called (after freeing the previous mem) whenever the image size changes.
 */
int alloc_man_mem(man_calc_struct* m, int width, int height)
{
	int n;

	m->iter_data_line_size = width + 2;
	m->image_size = width * height; // new image size

	// Because the fast algorithm checks offsets from the current pixel 
	// location, iter_data needs dummy lines to accomodate off-screen checks. 
	// Needs one line at y = -1, and 6 at y = ysize. Also needs two dummy 
	// pixels at the end of each line.

	// Need separate pointer to be able to free later

	m->iter_data_start = (unsigned int*)malloc
	(n = m->iter_data_line_size * (height + 7) * sizeof(m->iter_data_start[0]));
	
	if (m->iter_data_start != NULL)
		memset(m->iter_data_start, 0, n);

	// create dummy lines at y = -1 for fast alg
	m->iter_data = m->iter_data_start + m->iter_data_line_size;

	// Create a corresponding array for the magnitudes. Don't really need the 
	// dummy lines but this allows using a fixed offset from iter_data
	m->mag_data = (float*)malloc
	(m->iter_data_line_size * (height + 7) * sizeof(m->mag_data[0]));

	m->mag_data_offs = (int)((char*)m->mag_data - (char*)m->iter_data);

	// These two need 4 extra dummy values
	m->img_re = (double*)malloc((width + 4) * sizeof(m->img_re[0]));
	m->img_im = (double*)malloc((height + 4) * sizeof(m->img_im[0]));

	// Buffer for PNG save (not needed for main calculation). 4 bytes per pixel
	if (m->flags & FLAG_IS_SAVE)
	{
		m->png_buffer = (unsigned char*)malloc
		((width << 2) * height * sizeof(unsigned char));
		
		if (m->png_buffer == NULL)
			return 0;
	}

	if (m->iter_data_start == NULL || 
		m->mag_data == NULL || 
		m->img_re == NULL || 
		m->img_im == NULL)
		return 0;
	return 1;
}

/* Free all memory allocated above */
void free_man_mem(man_calc_struct* m)
{
	if (m->iter_data_start != NULL)
	{
		free(m->iter_data_start);
		free(m->mag_data);
		free(m->img_re);
		free(m->img_im);
		if (m->png_buffer != NULL)
			free(m->png_buffer);
	}
}

/**
 * Rename this; now does a lot more than create a bitmap
 *
 * Currently called every time the window is resized by even one pixel. Really 
 * should only be called when user is done with the resizing movement.
 *
 * Only called prior do doing main calculation (not for save)
 */
int create_bitmap(int width, int height)
{
	BITMAPINFO bmi;
	BITMAPINFOHEADER* h;
	int i, j, err;
	int bmihsize = sizeof(BITMAPINFOHEADER);
	static int prev_width = 0;
	static int prev_height = 0;
	man_calc_struct* m = &main_man_calc_struct;

	// If size didn't change, return
	if (prev_width == width && prev_height == height)
		return 0;

	// free any existing arrays/bitmaps
	if (m->iter_data_start != NULL)
		for (i = 0; i < 4; i++)
			DeleteObject(quad[i].handle);
	free_man_mem(m);

	memset(&bmi, 0, bmihsize);
	h = &bmi.bmiHeader;
	h->biSize = bmihsize;
	h->biWidth = width;
	h->biHeight = -height;
	h->biPlanes = 1;
	h->biBitCount = 32;
	h->biCompression = BI_RGB;

	// Create the 4 quadrant bitmaps
	err = 0;
	for (i = 0; i < 4; i++)
	{
		quad[i].handle = CreateDIBSection
		(NULL, &bmi, DIB_RGB_COLORS, (void**)&quad[i].bitmap_data, NULL, 0);

		if (!quad[i].handle)
			err = 1;
	}

	if (!hscreen_dc || err || 
		!alloc_man_mem(&main_man_calc_struct, width, height))
	{
		MessageBox(NULL, "Error allocating storage arrays.", 
			NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
		// really should exit cleanly here... subsequent code will crash
	}

	// update re/im from any pan offsets and reset offsets resized; need to 
	// recalculate entire image
	update_re_im(m, m->pan_xoffs, m->pan_yoffs);
	status |= STAT_NEED_RECALC;
	
	prev_width = width;
	prev_height = height;
	
	// set smaller dimension
	m->min_dimension = (width < height) ? width : height;

	// Precalculate pointer offsets of neighboring pixels for the fast "wave" 
	// algorithm (these only change when image width changes). Not needed for 
	// save
	for (i = 1; i < 7; i++)
		for (j = 0; j < 4; j++)
			wave_ptr_offs[i][j] = 
			wave_yoffs[i][j] * m->iter_data_line_size + wave_xoffs[i][j];

	reset_quadrants();	// reset to recalculate all
	reset_fps_values();	// reset frames/sec timing values
	reset_pan_state();	// reset pan filters and movement state

	return 1;
}

void init_combo_box(HWND hwnd, int dlg_item, char** strs, int n, 
	int default_selection)
{
	int i;
	// Initialize all the strings
	for (i = 0; i < n; i++)
		SendDlgItemMessage(hwnd, dlg_item, CB_ADDSTRING, 0, (LPARAM)strs[i]);

	// Set the default selection
	SendDlgItemMessage(hwnd, dlg_item, CB_SETCURSEL, default_selection, 0);
}

/* Find a string in a strings array and return its index, or -1 if not found */
int get_string_index(char* str, char** strs, int num_strs)
{
	int i;
	for (i = 0; i < num_strs; i++)
		if (!strcmp(str, strs[i]))
			return i;
	return -1;
}

/* Returns 1 if palette is one of the builtins, else 0 */
int get_builtin_palette(void)
{
	int tmp;
	char str[256];
	man_calc_struct* m = &main_man_calc_struct;

	GetDlgItemText(hwnd_dialog, IDC_PALETTE, str, sizeof(str));

	if ((tmp = get_string_index(str, palette_strs, NUM_ELEM(palette_strs))) >= 0)
	{
		m->palette = tmp;
		return 1;
	}
	return 0;
}

/**
 * Try loading a user palette from a file. Always reload so user can edit files
 * on the fly. If file is missing or bad, palette doesn't change.
 */
void get_user_palette(void)
{
	FILE* fp;
	unsigned int tmp, bmp_flag;
	man_calc_struct* m = &main_man_calc_struct;

	GetDlgItemText(hwnd_dialog, IDC_PALETTE, palette_file, sizeof(palette_file));

	// See whether it's a BMP or text palette file; load using corresponding 
	// function. This is safe because files have to have 3-char extensions to 
	// make it into the dropdown list.
	bmp_flag = !_strnicmp(palette_file + strlen(palette_file) - 3, "bmp", 3);

	if ((fp = open_file(palette_file, "", bmp_flag)) != NULL)
	{
		tmp = bmp_flag ? load_palette_from_bmp(fp, m) : load_palette(fp, m);

		if (tmp) // load_palette* assigns a nonzero number to the palette
			m->palette = tmp;	// if valid, which can then be used with 
								// apply_palette.
		else
		{
			MessageBox(NULL, bmp_flag
				? "Unsupported file format. "
				  "Please supply an uncompressed 24-bit bitmap."
				: "Unrecognized file format.",
				NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
		}
		fclose(fp);
	}
}

int get_rendering_alg(void)
{
	char str[256];

	GetDlgItemText(hwnd_dialog, IDC_RENDERING, str, sizeof(str));
	return get_string_index(str, rendering_strs, NUM_ELEM(rendering_strs));
}

int get_precision(void)
{
	char str[256];

	GetDlgItemText(hwnd_dialog, IDC_PRECISION, str, sizeof(str));
	return get_string_index(str, precision_strs, NUM_ELEM(precision_strs));
}

int get_alg(void)
{
	char str[256];

	GetDlgItemText(hwnd_dialog, IDC_ALGORITHM, str, sizeof(str));
	return get_string_index(str, alg_strs, NUM_ELEM(alg_strs));
}

void get_num_threads(void)
{
	char str[256];

	GetDlgItemText(hwnd_dialog, IDC_THREADS, str, sizeof(str));
	num_threads_ind = get_string_index
		(str, num_threads_strs, NUM_ELEM(num_threads_strs));
	
	num_threads = 1 << num_threads_ind;
}

/* Increase, decrease, or just clip the max iterations */
void update_iters(int up, int down)
{
	man_calc_struct* m = &main_man_calc_struct;

	if (up)
		m->max_iters <<= 1;
	if (down)
		m->max_iters >>= 1;
	if (m->max_iters < MIN_ITERS)
		m->max_iters = MIN_ITERS;
	if (m->max_iters > MAX_ITERS)
		m->max_iters = MAX_ITERS;

	SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
}

/**
 * Reset mandelbrot parameters to the home image
 */
void set_home_image(void)
{
	man_calc_struct* m = &main_man_calc_struct;

	m->pan_xoffs = 0;	// Reset any pan offsets
	m->pan_yoffs = 0;
	m->re = HOME_RE;
	m->im = HOME_IM;
	m->mag = HOME_MAG;
	m->max_iters = HOME_MAX_ITERS;	// Better to reset the max iters here. 
									// Don't want large #
	update_iters(0, 0);	// from previous image
}

/**
 * Get all fields from the dialog
 */
void get_dialog_fields(void)
{
	man_calc_struct* m = &main_man_calc_struct;

	m->max_iters = GetDlgItemInt(hwnd_dialog, IDC_ITERS, NULL, FALSE);
	update_iters(0, 0); // clip

	// Doing these only in the dialog box handler mean they don't take effect 
	// until user actually clicks them, whereas getting them here returns the 
	// realtime values.

	// It's useful to be able to set these on the fly
	m->alg = get_alg();
	m->precision = get_precision();
	m->rendering_alg = get_rendering_alg();

	// don't allow on-the-fly change to this (unimplemented)
	if (m->precision == PRECISION_EXTENDED)
		m->precision = PRECISION_DOUBLE;

	get_builtin_palette();

	// This not so much
	// get_num_threads();
}

/* Update the slider position and get the new value. Clips any bad values. */
int set_slider_pos(int dlg_item, int pos)
{
	SendDlgItemMessage(hwnd_dialog, dlg_item, TBM_SETPOS, TRUE, (LONG)pos);
	return (int)SendDlgItemMessage(hwnd_dialog, dlg_item, TBM_GETPOS, 0, 0);
}

void setup_sliders(void)
{
	cfg_default.pan_rate.c_val = set_slider_pos
		(IDC_PAN_RATE, cfg_default.pan_rate.c_val);
	
	cfg_default.zoom_rate.c_val = set_slider_pos
		(IDC_ZOOM_RATE, cfg_default.zoom_rate.c_val);
}

/**
 * Print status line at bottom of dialog with current precision, logfile image 
 * position, and logfile image total. Also print precision loss indicator if 
 * applicable.
 */
void print_status_line(int calc)
{
	char s[128];
	man_calc_struct* m = &main_man_calc_struct;

	// if currently saving, keep saving status in first part of line
	if (!(status & STAT_DOING_SAVE))
	{
		sprintf_s(s, sizeof(s), "%s%s", calc ? "Calculating..." : "Ready ",
			calc ? "" : precision_loss ? "[Prec Loss]" : "");

		SetWindowText(hwnd_status, s);
	}

	sprintf_s(s, sizeof(s), "%d/%d  %c", log_pos + 1, log_count,
		m->precision == PRECISION_SINGLE 
		? 'S' : m->precision == PRECISION_DOUBLE ? 'D' : 'E');
	SetWindowText(hwnd_status2, s);
}

/**
 * If the palette has been modified in some way (inverted, locked, colors 
 * changed with xor, etc) print a '*' before the Palette text.
 */
void print_palette_status(void)
{
	char s[32];
	man_calc_struct* m = &main_man_calc_struct;

	sprintf_s(s, sizeof(s), "%c Palette", 
		(status & STAT_PALETTE_LOCKED) || m->pal_xor ? '*' : ' ');
	SetWindowText(GetDlgItem(hwnd_dialog, IDC_PAL_TEXT), s);
}

void not_implemented_yet(void)
{
	MessageBox(NULL, "This feature is not implemented yet.", 
		NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
}

void unsupported_alg_prec(void)
{
	MessageBox(NULL, 
		"Your (obsolete) CPU cannot run this algorithm/precision combination.\n"
		"Using C algorithm.", NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
}

int unrecommended_alg(void)
{
	return MessageBox(NULL, 
		"Using the Fast algorithm with Normalized rendering may\n"
		"cause image artifacts. Switch to the Exact algorithm?",
		"Warning", MB_YESNO | MB_ICONWARNING | MB_TASKMODAL);
}

/**
 * Print ! before algorithm to indicate a warning if using normalized rendering
 * and fast alg.
 */
void set_alg_warning(void)
{
	char s[32];
	man_calc_struct* m = &main_man_calc_struct;

	sprintf_s(s, sizeof(s), "  Algorithm");
	if (!(m->alg & ALG_EXACT) && (m->rendering_alg == RALG_NORMALIZED))
		s[0] = '!';
	SetWindowText(GetDlgItem(hwnd_dialog, IDC_ALGORITHM_TEXT), s);
}

/**
 * See if the CPU supports the selected algorithm/precision combination, and 
 * whether the alg is suitable for normalized rendering mode.
 */
void check_alg(HWND hwnd)
{
	man_calc_struct* m = &main_man_calc_struct;

	m->precision = get_precision();
	m->alg = get_alg();
	m->rendering_alg = get_rendering_alg();

	if (m->precision == PRECISION_EXTENDED)
	{
		not_implemented_yet();
		SendDlgItemMessage(hwnd, IDC_PRECISION, CB_SETCURSEL, 
			PRECISION_DOUBLE, 0);
	}
	else if (m->precision == PRECISION_DOUBLE)
	{
		// If CPU doesn't support SSE2, can only run C version
		if (sse_support < 2 && !(m->alg & ALG_C))
		{
			unsupported_alg_prec();
			SendDlgItemMessage(hwnd, IDC_ALGORITHM, CB_SETCURSEL, 
				ALG_FAST_C, 0);
		}
	}
	else
	{
		// If CPU doesn't support SSE, can only run C version
		if (!sse_support && !(m->alg & ALG_C))
		{
			unsupported_alg_prec();
			SendDlgItemMessage(hwnd, IDC_ALGORITHM, CB_SETCURSEL, 
				ALG_FAST_C, 0);
		}
	}

	// Give user a choice to switch to exact alg if using normalized rendering

	if (!(m->alg & ALG_EXACT) && (m->rendering_alg == RALG_NORMALIZED))
		if (unrecommended_alg() == IDYES)
		{
			SendDlgItemMessage
			(hwnd, IDC_ALGORITHM, CB_SETCURSEL, m->alg |= ALG_EXACT, 0);
			
			// need to recalc if switching to exact
			status |= STAT_RECALC_FOR_PALETTE;
		}
	set_alg_warning();
}

/**
 * Calculate all or a portion of the set. If recalc_all is nonzero, recalculate 
 * entire image, update image info/status line and set a wait cursor during 
 * calculation. Otherwise only recalculate the update rectangles and don't 
 * update info.
 *
 * sdmc
 */
void do_man_calculate(int recalc_all)
{
	HCURSOR cursor = NULL;
	man_calc_struct* m = &main_man_calc_struct;

	// make max iters even (required by optimized algorithm)
	m->max_iters &= ~1;
	
	// need to recalculate all if max iters changed
	if (m->max_iters != m->max_iters_last)
		status |= STAT_NEED_RECALC;
	if (status & STAT_NEED_RECALC)	// if need recalculation,
		recalc_all = 1;	// force info update and wait cursor
	if (recalc_all)
	{
		// could be recalculating for palette- renormalize
		update_re_im(m, m->pan_xoffs, m->pan_yoffs);
		
		reset_quadrants();	// reset quadrants to UL only
		if (!do_rtzoom)		// If not realtime zooming
		{
			// print status line and reset fps timing values
			print_status_line(1);
			reset_fps_values();
			cursor = GetCursor();	// save current cursor and set wait cursor
			SetCursor(wait_cursor);
		}
		// do this always, so we can change max_iters while zooming, etc.
		get_dialog_fields();
		
		// no longer need to recalc for palette
		status &= ~STAT_RECALC_FOR_PALETTE;
	}

	man_calculate_quadrants();
	
	// last iters actually calculated, for palette code
	m->max_iters_last = m->max_iters;

	InvalidateRect(hwnd_main, NULL, 0); // cause repaint with image data
	UpdateWindow(hwnd_main);

	// Don't update this stuff if realtime zooming. Will be done at intervals
	if (recalc_all && !do_rtzoom)
	{
		SetWindowText(hwnd_info, get_image_info(1));	// update time, GFlops
		print_status_line(0);
		SetCursor(cursor);	// restore old cursor
	}
}

/**
 * Need to calculate all this stuff to accomodate strange window borders, large 
 * fonts and such.
 */
void get_system_metrics(void)
{
	x_border = 2 * GetSystemMetrics(SM_CXSIZEFRAME);
	y_border = 2 * (y_thinborder = GetSystemMetrics(SM_CYSIZEFRAME)) + 
		GetSystemMetrics(SM_CYCAPTION);

	// X/Y reversed on these
	x_dialog_border = 2 * GetSystemMetrics(SM_CXFIXEDFRAME);
	y_dialog_border = 2 * GetSystemMetrics(SM_CYFIXEDFRAME) + 
		GetSystemMetrics(SM_CYSMCAPTION);

	// Get font size. Assumes x/y pixels per inch are the same. Returns 96 for 
	// normal, 120 for large.
	lpix_per_inch = GetDeviceCaps(GetDC(NULL), LOGPIXELSX);
}

/**
 * Update the control dialog. If move is nozero, moves it to the edge of the 
 * mandelbrot window. If hide is nonzero, hides it, otherwise shows it.
 */
void update_dialog(int hide, int move)
{
	HWND hwnd_desktop;
	RECT rc_owner, rc_dialog, rc_desktop;
	int xpos, ypos, overhang;

	// this can be called before the main window gets created
	if (hwnd_main == NULL)
		return;

	GetWindowRect(hwnd_dialog, &rc_dialog);

	if (move)
	{
		hwnd_desktop = GetDesktopWindow();
		GetWindowRect(hwnd_main, &rc_owner);
		GetWindowRect(hwnd_desktop, &rc_desktop);

		xpos = rc_owner.right;
		ypos = rc_owner.top;

		// Clip so dialog doesn't go off the right end of the desktop. Also 
		// move down so minimize/maximize buttons are still visible. Dialog 
		// right - left = width
		overhang = xpos + (rc_dialog.right - rc_dialog.left) - rc_desktop.right;

		if (overhang > 0)
		{
			xpos -= overhang;
			ypos += y_border - y_thinborder;
		}
	}
	else // keep at current position
	{
		xpos = rc_dialog.left;
		ypos = rc_dialog.top;
	}

	SetWindowPos(hwnd_dialog, HWND_TOP, xpos, ypos, 0, 0, SWP_NOSIZE |
		(hide ? SWP_HIDEWINDOW : SWP_SHOWWINDOW));
}

/* Enter and exit fullscreen mode. */
void toggle_fullscreen(void)
{
	if ((status ^= STAT_FULLSCREEN) & STAT_FULLSCREEN)
	{
		// Hide dialog on entry to fullscreen mode, if so configured
		if (!(cfg_default.options.c_val & OPT_DIALOG_IN_FULLSCREEN))
			status |= STAT_DIALOG_HIDDEN;

		SetWindowLongPtr(hwnd_main, GWL_STYLE, WS_POPUP | WS_VISIBLE);
		SetWindowPos(hwnd_main, NULL, 0, 0, GetSystemMetrics(SM_CXSCREEN),
			GetSystemMetrics(SM_CYSCREEN), SWP_DRAWFRAME | SWP_NOZORDER);
	}
	else // exit
	{
		// always show dialog when coming out of fullscreen
		status &= ~STAT_DIALOG_HIDDEN;

		SetWindowLongPtr(hwnd_main, GWL_STYLE, 
			WS_OVERLAPPEDWINDOW | WS_VISIBLE);
		SetWindowPos(hwnd_main, NULL, main_rect.left, main_rect.top,
			main_rect.right - main_rect.left, 
			main_rect.bottom - main_rect.top,
			SWP_NOZORDER);
	}
	UpdateWindow(hwnd_main);

	// recalc after resize, if enabled
	if (cfg_default.options.c_val & OPT_RECALC_ON_RESIZE)
		status |= STAT_RECALC_IMMEDIATELY;
}

unsigned int CALLBACK show_help(LPVOID param)
{
	// Using MB_TOPMOST gives an old window style... need to use SYSTEMMODAL, 
	// but then we get an annoying icon on the title bar.
	MessageBox(NULL, help_text, "QuickMAN Help", 
		MB_OK /*| MB_ICONINFORMATION*/ | MB_SYSTEMMODAL);
	status &= ~STAT_HELP_SHOWING;
	return 0;
}

/**
 * Resize window if necessary. Don't allow window size changes in fullscreen 
 * mode, except by restore command (assume if user went to fullscreen, he wants 
 * to stay there).
 *
 * v1.07
 */
void resize_window(void)
{
	if (!(status & STAT_FULLSCREEN))
	{
		// xsize 0: fullscreen (other values < MIN_SIZE can be used later)
		if (cfg_default.xsize.c_val < MIN_SIZE)
			toggle_fullscreen();
		else if (cfg_default.ysize.c_val >= MIN_SIZE && // ignore any restoring ysize = 0
				(cfg_default.xsize.c_val != prev_xsize || 
				cfg_default.ysize.c_val != prev_ysize))
		{
			// Without this, you can change the size of a maximized window, 
			// which puts it into a bad state (can't resize)
			ShowWindow(hwnd_main, SW_RESTORE);
			SetWindowPos(hwnd_main, HWND_TOP, 0, 0, 
				cfg_default.xsize.c_val + x_border,
				cfg_default.ysize.c_val + y_border, SWP_NOMOVE | SWP_NOCOPYBITS);

			// not sure why this is necessary. Sizes get out of sync?
			UpdateWindow(hwnd_main);

			prev_xsize = cfg_default.xsize.c_val;	// Save previous size
			prev_ysize = cfg_default.ysize.c_val;	// Can't do this in 
													// create_bitmap
		}
	}
	else if (cfg_default.ysize.c_val < MIN_SIZE)
		toggle_fullscreen(); // If in fullscreen mode: only allow restore
}

/**
 * Save the image in one or more formats: PNG image, pixel iteration counts 
 * (32 bit unsigneds), and/or pixel final magnitudes (32-bit floats). Only the 
 * PNG save is currently implemented; selection checkboxes are grayed out.
 *
 * This runs in its own thread, so it can be happening in the background during 
 * normal browsing.
 *
 * Currently only does one row at a time. A more complex chunk-based version is 
 * only about 3% faster and less friendly as a background task.
 */
unsigned int CALLBACK do_save(LPVOID param)
{
	int i, j, n, save_xsize, save_ysize;
	unsigned char* ptr3, * ptr4, c[256];
	FILE* fp;

	TIME_UNIT start_time, t;
	man_calc_struct* m, * s;

	m = &main_man_calc_struct;
	s = &save_man_calc_struct;

	// Get sizes for image save, and clip
	if ((save_xsize = GetDlgItemInt(hwnd_dialog, 
		IDC_SAVE_XSIZE, NULL, FALSE)) < MIN_SIZE)
		save_xsize = MIN_SIZE;
	if ((save_ysize = GetDlgItemInt(hwnd_dialog, 
		IDC_SAVE_YSIZE, NULL, FALSE)) < MIN_SIZE)
		save_ysize = MIN_SIZE;

	// Get filename and add .png extension if not already present
	GetDlgItemText(hwnd_dialog, IDC_SAVEFILE, savefile, sizeof(savefile));
	n = (int)strlen(savefile);
	if (n < 4 || _strnicmp(&savefile[n - 4], ".png", 4))
		strcat_s(savefile, sizeof(savefile), ".png");

	s->xsize = save_xsize;
	s->ysize = 1;

	// Copy relevant image parameters to the save calculation structure from 
	// the main structure (whose image is currently displayed).

	// Get saved image re/im from the main re/im + pan offsets
	s->re = m->re + get_re_im_offs(m, m->pan_xoffs);
	s->im = m->im - get_re_im_offs(m, m->pan_yoffs);

	// min_dimension is used to calc. coords- needs to reflect the actual ysize
	s->min_dimension = (save_xsize > save_ysize) ? save_ysize : save_xsize;
	s->mag = m->mag;
	s->max_iters = s->max_iters_last = m->max_iters;

	// Can get unexpected precision loss when the saved image is larger than 
	// the on-screen image.
	// Always use best precision to minimize occurrences
	s->precision = PRECISION_DOUBLE;	// m->precision
	s->alg = m->alg | ALG_EXACT;		// exact will be faster for 1-pixel 
										// high rows. Want for best quality 
										// anyway.
	s->palette = m->palette;
	s->prev_pal = 0xFFFFFFFF;			// always recalc. pal lookup table 
										// before starting
	s->pal_xor = m->pal_xor;
	s->max_iters_color = m->max_iters_color;
	s->rendering_alg = m->rendering_alg;
	s->flags |= FLAG_CALC_RE_ARRAY;		// tell man_calculate to calculate the 
										// real array initially

	// Make sure all image data above is already captured before possibly 
	// popping a message box. Because this func is in a separate thread, the 
	// image can be modified (panned, etc) while the box is open.

	// If file already exists, ask user to confirm overwite
	if (!fopen_s(&fp, savefile, "rb"))
	{
		fclose(fp);
		sprintf_s(c, sizeof(c), "%s already exists. Overwrite?", savefile);
		if (MessageBox(NULL, c, "Warning", 
			MB_YESNO | MB_ICONWARNING | MB_TASKMODAL) != IDYES)
		{
			status &= ~STAT_DOING_SAVE;
			return 0;
		}
	}

	if (!png_save_start(savefile, save_xsize, save_ysize))
	{
		status &= ~STAT_DOING_SAVE;
		return 0;
	}

	free_man_mem(s); // free any existing arrays and alloc new arrays
	alloc_man_mem(s, save_xsize, 1);

	start_time = t = get_timer();

	s->pan_yoffs = -((save_ysize - 1) >> 1); // pan_yoffs of image top
	for (i = 0; i < save_ysize; i++)
	{
		man_calculate(s, 0, save_xsize - 1, 0, 0); // iterate the row
		
		// don't have to recalculate real array on subsequent rows
		s->flags &= ~FLAG_CALC_RE_ARRAY;

		// Palette-map the iteration counts to RGB data in png_buffer. 
		// Magnitudes are also available here.
		apply_palette(s, (unsigned int*)s->png_buffer, s->iter_data, 
			save_xsize, 1, num_threads);

		// Convert the 4 bytes-per-pixel data in png_buffer to 3 bpp, as 
		// required by pnglib, and write rows
		ptr3 = ptr4 = s->png_buffer;
		for (j = 0; j < save_xsize; j++)
		{
			*((unsigned int*)ptr3) = *((unsigned int*)ptr4);
			ptr3 += 3;
			ptr4 += 4;
		}
		// write the row
		if (!png_save_write_row(s->png_buffer))
			break;

		s->pan_yoffs++; // go to next row

		// print progress indicator every 0.5s
		if (get_seconds_elapsed(t) > 0.5)
		{
			sprintf_s(c, sizeof(c), "Saving... (%3.1f%%)", 
				100.0 * (double)i / (double)save_ysize);
			SetWindowText(hwnd_status, c);
			t = get_timer();
		}
	}

	png_save_end();
	sprintf_s(c, sizeof(c), "Saved in %.1fs", get_seconds_elapsed(start_time));
	SetWindowText(hwnd_status, c);

	status &= ~STAT_DOING_SAVE;
	return 1;
}

/**
 * Handler for the control dialog box. How do we make this stop blocking the 
 * rest of the code when the user is moving it?
 *
 * sdiag
 */
INT_PTR CALLBACK 
man_dialog_proc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	static int adj_iters_prev = 0, ignore_next_change = 0;
	static int new_file_entered = 0, new_file_selected = 0;
	int dialog_w, dialog_h, text_w, text_h, frame_w, frame_h, last_ypixel;
	int done, max_tries;
	int tab_spacing = 26; // dialog box units
	RECT rc_tframe, rc_status, rc_dialog;
	log_entry* e;
	man_calc_struct* m = &main_man_calc_struct;

	// Original dialog box coords and sizes, for adjustment
#define ORIG_YBORDER		34
#define ORIG_DIALOG_HEIGHT	(700 + ORIG_YBORDER)
#define ORIG_DIALOG_WIDTH	173
#define ORIG_STATUS_X		8
#define YBORDER_ADJUSTMENT	(y_border - ORIG_YBORDER)
#define ORIG_LPIX			96	// original logical pixels per inch
#define TFRAME_MIN_HEIGHT	1	// thumbnail frame is now an invisible dummy 
								// item used for sizing

	switch (uMsg)
	{
	case WM_INITDIALOG:

		// Initialize all dialog box fields
		SetDlgItemInt(hwnd, IDC_ITERS, m->max_iters, FALSE);
		SendDlgItemMessage(hwnd, IDC_PAN_RATE, TBM_SETRANGE, TRUE, 
			MAKELONG(0, MAX_PAN_RATE));
		SendDlgItemMessage(hwnd, IDC_ZOOM_RATE, TBM_SETRANGE, TRUE, 
			MAKELONG(0, MAX_ZOOM_RATE));
		init_combo_box(hwnd, IDC_PRECISION, precision_strs, 
			NUM_ELEM(precision_strs), m->precision);
		init_combo_box(hwnd, IDC_PALETTE, palette_strs, 
			NUM_ELEM(palette_strs), m->palette);
		init_combo_box(hwnd, IDC_RENDERING, rendering_strs, 
			NUM_ELEM(rendering_strs), m->rendering_alg);
		init_combo_box(hwnd, IDC_ALGORITHM, alg_strs, 
			NUM_ELEM(alg_strs), m->alg);
		init_combo_box(hwnd, IDC_THREADS, num_threads_strs, 
			MAX_THREADS_IND + 1, num_threads_ind);
		init_combo_box(hwnd, IDC_LOGFILE, file_strs, 
			NUM_ELEM(file_strs), 0);

		// default save filename- later have this scan for next available
		SetWindowText(GetDlgItem(hwnd, IDC_SAVEFILE), savefile);

		// Set tab stops for the INFO window
		SendDlgItemMessage(hwnd, IDC_INFO, EM_SETTABSTOPS, 1, 
			(LPARAM)&tab_spacing);

		// If any of these are null, the system has major problems
		hwnd_iters = GetDlgItem(hwnd, IDC_ITERS);
		hwnd_info = GetDlgItem(hwnd, IDC_INFO);
		hwnd_status = GetDlgItem(hwnd, IDC_STATUS);
		hwnd_status2 = GetDlgItem(hwnd, IDC_STATUS2);
		// v1.10: now an invisible dummy item
		hwnd_thumbnail_frame = GetDlgItem(hwnd, IDC_THUMBNAIL_FRAME);

		// PNG save and Preserve aspect ratio buttons checked by default
		SendDlgItemMessage(hwnd, IDC_PNG, BM_SETCHECK, BST_CHECKED, 0);
		SendDlgItemMessage(hwnd, IDC_ASPECT, BM_SETCHECK, BST_CHECKED, 0);

		// Adjust dialog box to compensate for different border and font sizes. 
		// What a pain...

		// To simulate Vista large fonts on XP: set font size to large, set 
		// lpix_per_inch to 120, and change "MS Shell Dlg" font size to 10 in 
		// quickman.rc. Seems almost exact.
		// lpix_per_inch = 120;

		// The seems to be no documented function to help figure out how a 
		// dialog box width should scale with font size. The calc. below 
		// usually turns out too wide. Will be corrected in an iterative 
		// process below...

		dialog_w = (lpix_per_inch * (ORIG_DIALOG_WIDTH - x_dialog_border)) / 
			ORIG_LPIX + x_dialog_border;

		// Initially, see if things fit in a dialog that's the same height as 
		// the main window
		dialog_h = ORIG_DIALOG_HEIGHT + YBORDER_ADJUSTMENT;

		// With pathologically large fonts, could fail to make a decent box- 
		// just give up after a certain number of tries
		max_tries = 3;

		do
		{
			SetWindowPos(hwnd, HWND_TOP, 0, 0, dialog_w, dialog_h, 
				SWP_NOMOVE | SWP_HIDEWINDOW);
			GetWindowRect(hwnd, &rc_dialog); // get adjusted dialog rect

			// Dialog item positions are relative to dialog box coords. y = 0 
			// is the topmost pixel in the dialog (inside the border). 
			//// Last_ypixel is the bottommost pixel in the dialog. These seem
			// to be the right calculations...

			last_ypixel = rc_dialog.bottom - rc_dialog.top - y_dialog_border - 1;

			GetWindowRect(hwnd_status, &rc_status);		// get status line rect
			
			// height of status line text
			text_h = rc_status.bottom - rc_status.top;

#define TEXT_TO_BOTTOM_SPACE  8 // some fine-tuning constants
#define TEXT_TO_FRAME_SPACE   8
#define FRAME_TO_DIALOG_SPACE 8 // distance to dialog edges from each vertical 
								// side of frame

			text_h += TEXT_TO_BOTTOM_SPACE;

			// Adjust status1/status2 ypos so they're just above the bottom of 
			// the dialog
			SetWindowPos(hwnd_status, HWND_TOP, ORIG_STATUS_X, 
				last_ypixel - text_h, 0, 0, SWP_NOSIZE);

			GetWindowRect(hwnd_status2, &rc_status);	// get status2 line rect
			
			// not sure why this width is off by a few pixels
			text_w = rc_status.right - rc_status.left + 4;

			// Get thumbnail frame rectangle
			GetWindowRect(hwnd_thumbnail_frame, &rc_tframe);

			// Adjust status2 xpos so its right lines up with the thumbnail 
			// frame right
			SetWindowPos(hwnd_status2, HWND_TOP, 
				rc_tframe.right - rc_dialog.left - text_w,
				last_ypixel - text_h, 0, 0, SWP_NOSIZE);

			// Adjust height of thumbnail frame so its bottom is just above the 
			// status line top
			SetWindowPos(hwnd_thumbnail_frame, HWND_TOP, 0, 0, 
				frame_w = rc_tframe.right - rc_tframe.left,
				frame_h = rc_dialog.bottom - rc_tframe.top - text_h - TEXT_TO_FRAME_SPACE,
				SWP_NOMOVE);

			done = 1;

			// If dialog is too wide or too narrow, adjust to frame width and 
			// start over
			if (dialog_w != (frame_w += 2 * (lpix_per_inch * FRAME_TO_DIALOG_SPACE) / 
				ORIG_LPIX + x_dialog_border))
			{
				dialog_w = frame_w;
				done = 0;
			}

			// If frame height is too small, increase dialog height and start 
			// over
			if (frame_h < TFRAME_MIN_HEIGHT)
			{
				dialog_h += TFRAME_MIN_HEIGHT - frame_h;
				done = 0;
			}
		} while (!done && --max_tries);
		return FALSE;

	case WM_VSCROLL: // Update iterations from spin control
		if ((HWND)lParam == GetDlgItem(hwnd, IDC_ADJUST_ITERS))
		{
			// Get value user may have edited
			m->max_iters = GetDlgItemInt(hwnd, IDC_ITERS, NULL, FALSE);
			switch (LOWORD(wParam))
			{
			case SB_THUMBPOSITION:
				// No change (0) or decreasing value means up, increasing value 
				// means down
				if (HIWORD(wParam) > adj_iters_prev)
					update_iters(0, 1);
				else
					update_iters(1, 0);
				adj_iters_prev = HIWORD(wParam);
				break;
			}
		}
		return TRUE;

	case WM_HSCROLL:	// Slider values are being updated. Moving these around 
						// during zooming annoyingly creates about a 1% 
						// slowdown, even with no processing.

		if ((HWND)lParam == GetDlgItem(hwnd, IDC_PAN_RATE))
			cfg_default.pan_rate.c_val = (int)SendDlgItemMessage
				(hwnd, IDC_PAN_RATE, TBM_GETPOS, 0, 0);

		if ((HWND)lParam == GetDlgItem(hwnd, IDC_ZOOM_RATE))
			cfg_default.zoom_rate.c_val = (int)SendDlgItemMessage
				(hwnd, IDC_ZOOM_RATE, TBM_GETPOS, 0, 0);
		
		reset_fps_values(); // reset frames/sec timing values when pan or zoom 
							// rate changes
		return TRUE;

	case WM_COMMAND:

		switch (LOWORD(wParam))
		{
		case IDC_LOGFILE:
			if (HIWORD(wParam) == CBN_EDITCHANGE)
			{
				// Mark that user entered a new filename: will be added to box 
				// when a log func is used.
				new_file_entered = 1;	
				new_file_selected = 0;
			}
			// No longer really need to read logfile immediately on selection 
			// change, with quickman.cfg now available for default settings. 
			// Old Microsoft code was dangerous (no limit to string copy size).
			if (HIWORD(wParam) == CBN_SELCHANGE)
			{
				new_file_entered = 0;
				new_file_selected = 1;
			}
			return TRUE;

		case IDC_PALETTE:
		case IDC_RENDERING:
			if (HIWORD(wParam) == CBN_SELCHANGE)
			{
				// If palette is not one of the builtins, try loading from file
				if (!get_builtin_palette())
					get_user_palette();

				// Give warning if fast alg and normalized rendering; allow 
				// user to switch to exact first
				if (LOWORD(wParam) == IDC_RENDERING)
					check_alg(hwnd);

				// Recalculate all first if we need to
				if ((status & STAT_RECALC_FOR_PALETTE) || 
					(m->max_iters != m->max_iters_last))
				{
					// update re/im from any pan offsets and reset offsets
					update_re_im(m, m->pan_xoffs, m->pan_yoffs);
					do_man_calculate(1);
				}

				// Apply palette to the whole image (in UL quadrant here)
				apply_palette(m, quad[UL].bitmap_data, m->iter_data, 
					m->xsize, m->ysize, num_threads);
				InvalidateRect(hwnd_main, NULL, 0); // cause repaint
				UpdateWindow(hwnd_main);
			}
			return TRUE;

		case IDC_ALGORITHM:
		case IDC_PRECISION:
			if (HIWORD(wParam) == CBN_SELCHANGE)
				// make sure CPU supports alg/precision combination
				check_alg(hwnd);
			return TRUE;

		case IDC_THREADS:
			if (HIWORD(wParam) == CBN_SELCHANGE)
				get_num_threads();
			return TRUE;

		case IDC_SAVE_XSIZE:
		case IDC_SAVE_YSIZE:

			// Set save xsize or ysize automatically from the other if preserve 
			// aspect box is checked
			if (HIWORD(wParam) == EN_UPDATE)
			{
				if (IsDlgButtonChecked(hwnd, IDC_ASPECT))
				{
					double aspect;
					aspect = (double)m->xsize / (double)m->ysize;

					// Need the ignore_ logic to prevent infinite message loops
					if (LOWORD(wParam) == IDC_SAVE_XSIZE)
						if (ignore_next_change != 1)
						{
							ignore_next_change = 2;
							SetDlgItemInt(hwnd, IDC_SAVE_YSIZE, 
								(int)(0.5 + (double)GetDlgItemInt
								(hwnd, IDC_SAVE_XSIZE, NULL, FALSE) / aspect), 
								FALSE);
						}
						else
							ignore_next_change = 0;
					else
						if (ignore_next_change != 2)
						{
							ignore_next_change = 1;
							SetDlgItemInt(hwnd, IDC_SAVE_XSIZE, 
								(int)(0.5 + (double)GetDlgItemInt
								(hwnd, IDC_SAVE_YSIZE, NULL, FALSE) * aspect), 
								FALSE);
						}
						else
							ignore_next_change = 0;
				}
			}
			return TRUE;

		case ID_HOME:
			// Reset to base image coordinates; deliberate fallthru
			set_home_image();
			
			// only do this on home or log next/prev, not recalculation
			autoreset_settings(&cfg_default);

			resize_window();

		case ID_CALCULATE:

			do_rtzoom = prev_do_rtzoom = 0;	// stop any rt zoom in progress
			
			// update re/im from any pan offsets and reset offsets
			update_re_im(m, m->pan_xoffs, m->pan_yoffs);
			reset_pan_state();		// reset pan filters and movement state
			get_pan_steps(NULL, NULL, 0);	// reset any pan lock
			print_palette_status();
			do_man_calculate(1);	// calculate all pixels, update image info
			SetFocus(hwnd_main);	// allow arrow keys to work immediately for 
									// panning
			// set focus AFTER calculating, or you get an annoying blink
			return TRUE;

		case ID_LOG_IMAGE:
		case ID_LOG_PREV:
		case ID_LOG_NEXT:

			// Get the current filename
			GetDlgItemText(hwnd, IDC_LOGFILE, logfile, sizeof(logfile));

			// If it's new, read it. If it was a newly entered filename, add it 
			// to the list.
			if (new_file_entered || new_file_selected)
			{
				if (LOWORD(wParam) != ID_LOG_IMAGE)
				{
					log_read(logfile, "", 1);	// read logfile into array
					
					// for testing load-balancing alg
					reset_thread_load_counters();
				}
				if (new_file_entered)
					SendDlgItemMessage(hwnd, IDC_LOGFILE, CB_ADDSTRING, 0, 
						(LPARAM)logfile);
			}
			if (LOWORD(wParam) == ID_LOG_IMAGE)
			{
				// Need to update re/im with panning offsets, or logged 
				// coordinates will be wrong!
				// Added this line for v1.03 bug fix
				update_re_im(m, m->pan_xoffs, m->pan_yoffs);
				
				// reset pos if new file
				log_update(logfile, new_file_entered | new_file_selected);
				
				// update current/total number of log images
				print_status_line(0);
				
				// Ok if this stays till next calculation
				SetWindowText(hwnd_status, "Logged");
			}
			new_file_entered = 0;
			new_file_selected = 0;

			// If user wants a new log image...
			if ((LOWORD(wParam) == ID_LOG_NEXT || 
				LOWORD(wParam == ID_LOG_PREV)) && log_count)
			{
				// Autoreset any previous settings that need it
				autoreset_settings(&cfg_default);	
				
				// Clear any pan offsets (no need to update re/im as they are 
				// reset from logfile)
				m->pan_xoffs = 0;
				m->pan_yoffs = 0;
				
				// get next/prev image from logfile array
				if ((e = log_get(LOWORD(wParam) == ID_LOG_NEXT)) == NULL)
					return TRUE;

				// Update any new settings: 0 = don't copy to default_val
				copy_changed_settings(&cfg_default, &e->log_settings, 0);

				// Update sliders, info box, iters, and palette
				setup_sliders();
				update_iters(0, 0);
				UpdateWindow(hwnd_iters);

				do_rtzoom = prev_do_rtzoom = 0; // stop any rt zoom in progress
				reset_pan_state();	// reset pan filters and movement state
				
				// set pan lock (or turn it off if pan_key is 0)
				get_pan_steps(NULL, NULL, cfg_default.pan_key.c_val);

				// Change palette, max iters color, and inversion status if not 
				// locked
				if (!(status & STAT_PALETTE_LOCKED))
				{
					m->pal_xor = cfg_default.pal_xor.c_val;
					m->max_iters_color = cfg_default.max_iters_color.c_val;

					SendDlgItemMessage(hwnd, IDC_PALETTE, CB_SETCURSEL, 
						m->palette, 0);

					// if user palette, read it in
					if (m->palette >= num_builtin_palettes)
						get_user_palette();
				}

				print_palette_status();
				resize_window();

				// If not zooming, just calculate
				if (!cfg_default.zoom_in_out.c_val)
				{
					// calculate all pixels, update image info
					do_man_calculate(1);

					// already calculated: don't need to do again
					status &= ~STAT_RECALC_IMMEDIATELY;
					
					SetFocus(hwnd_main);
					return TRUE;
				}
				else
				{
					// Else fallthrough to do realtime zoom in (later need zoom 
					// out also)
				}
			}
			else
				return TRUE;

		case ID_ZOOM: // Do a realtime zoom in to the current image
			// update re/im from any pan offsets and reset offsets
			update_re_im(m, m->pan_xoffs, m->pan_yoffs);

			reset_fps_values();
			reset_thread_load_counters();
			zoom_start_time = get_timer();
			zoom_start_mag = m->mag;
			m->mag = MAG_MIN;
			do_rtzoom = RTZOOM_IN | RTZOOM_WITH_BUTTON;
			return TRUE;

		case ID_FULLSCREEN:
			toggle_fullscreen();
			return TRUE;

			// case ID_SLIDESHOW:
			// case ID_OPTIONS:
			   // not_implemented_yet();
			   // return TRUE;

		case ID_SAVE_IMAGE:
			if (!(status & STAT_DOING_SAVE))
			{
				// prevent re-entering function when already saving
				status |= STAT_DOING_SAVE;

				// With just WT_EXECUTEDEFAULT here (by accident), got strange 
				// behavior- sometimes wouldn't save (created a 0K file with no 
				// error indication or status update).
				// Hard to reproduce. With these two, seems ok. 
				// WT_EXECUTEINIOTHREAD necessary?
				QueueUserWorkItem(do_save, NULL, 
					WT_EXECUTELONGFUNCTION | WT_EXECUTEINIOTHREAD |
					(MAX_QUEUE_THREADS << 16));
			}
			return TRUE;

		case ID_HELP_BUTTON: // MS won't generate ID_HELP

			// Have to do this in a separate thread just to keep the messagebox 
			// from blocking my main window...

			// don't show help box if it's already showing
			if (!(status & STAT_HELP_SHOWING))
			{
				status |= STAT_HELP_SHOWING;
				QueueUserWorkItem(show_help, NULL, 
					WT_EXECUTELONGFUNCTION | (MAX_QUEUE_THREADS << 16));
			}
			return TRUE;

		default:
			return FALSE;
		}
		break;

	case WM_CLOSE:	// Never close or destroy this dialog
		break;
	case WM_DESTROY:
		break;
	}
	return FALSE;	// return FALSE if we didn't process the message 
					// (exceptions for INITDIALOG) otherwise TRUE
}

/**
 * Blits thin horizontal stripes from quadrant bitmaps (each stripe potentially
 * coming from 2 quadrants), so pixels are copied in upper-left to lower-right 
 * order.
 * Emulates blitting a single bitmap to avoid visual artifacts. Will be an 
 * exact emulation (with a lot of overhead) if STRIPE_THICKNESS is 1.
 */
void striped_blit(quadrant* ql, quadrant* qr, HDC hdc, HDC hscreen_dc)
{
	int src_yoffs, dest_yoffs, ysize, this_y, y_done;

	// Thickness of the stripes: the thinner the better, but thinner
	// stripes cause more overhead. 8 gives no measureable overhead on
	// the Athlon 4000+ system, but significant overhead on the Pentium
	// D 820 system. Probably hugely dependent on the video driver.

	// No artifacts visible with either 8 or 16 (except those that are present
	// with a full bitmap also- tearing, CPU cycle stealing by other 
	// applications, and frame rate aliasing with the screen refresh rate).

	//#define STRIPE_THICKNESS 16 // this now comes from global config setting

	// Return if no data in these quadrants
	if (!(ql->status & QSTAT_DO_BLIT) && !(qr->status & QSTAT_DO_BLIT))
		return;

	// The src_yoffs, dest_yoffs, and blit_ysize fields of the left and right
	// quadrants will always be the same.

	// Get from left quad if it has a blit rectangle, else get from right
	if (ql->status & QSTAT_DO_BLIT)
	{
		src_yoffs = ql->src_yoffs;
		dest_yoffs = ql->dest_yoffs;
		ysize = ql->blit_ysize;
	}
	else
	{
		src_yoffs = qr->src_yoffs;
		dest_yoffs = qr->dest_yoffs;
		ysize = qr->blit_ysize;
	}

	this_y = cfg_default.blit_stripe_thickness.c_val; // STRIPE_THICKNESS;
	y_done = 0;
	do
	{
		if (y_done + this_y > ysize)
			this_y = ysize - y_done;

		// Blit stripe left half from left quadrant
		if (ql->status & QSTAT_DO_BLIT)
		{
			SelectObject(hscreen_dc, ql->handle);
			BitBlt(hdc, ql->dest_xoffs, dest_yoffs, ql->blit_xsize, this_y,
				hscreen_dc, ql->src_xoffs, src_yoffs, SRCCOPY);
		}
		// Blit stripe right half from right quadrant
		if (qr->status & QSTAT_DO_BLIT)
		{
			SelectObject(hscreen_dc, qr->handle);
			BitBlt(hdc, qr->dest_xoffs, dest_yoffs, qr->blit_xsize, this_y,
				hscreen_dc, qr->src_xoffs, src_yoffs, SRCCOPY);
		}

		src_yoffs += this_y;
		dest_yoffs += this_y;
		y_done += this_y;
	} while (y_done != ysize);
}

/* Microsoft code for confining the mouse cursor to the main window. */
void confine_mouse_cursor(void)
{
	RECT rc;		// working rectangle
	POINT ptCUL;	// client upper left corner
	POINT ptCLR;	// client lower right corner

	// Retrieve the screen coordinates of the client area, and convert them 
	// into client coordinates.
	GetClientRect(hwnd_main, &rc);
	ptCUL.x = rc.left;
	ptCUL.y = rc.top;
	// Add one to the right and bottom sides, because the coordinates retrieved 
	// by GetClientRect do not.
	ptCLR.x = rc.right + 1;
	ptCLR.y = rc.bottom + 1;
	
	// include the far left and lowermost pixels.
	ClientToScreen(hwnd_main, &ptCUL);
	ClientToScreen(hwnd_main, &ptCLR);

	SetRect(&rc, 
		ptCUL.x, ptCUL.y,  // Copy the client coordinates of the client area
		ptCLR.x, ptCLR.y); // to the rcClient structure.

	SetCapture(hwnd_main);  // capture mouse input
	ClipCursor(&rc);        // confine the mouse cursor to the client area
}

/* The window function for the main window */
LRESULT CALLBACK 
MainWndProc(HWND hwnd, UINT nMsg, WPARAM wParam, LPARAM lParam)
{
	RECT rc;
	HDC hdc;
	PAINTSTRUCT ps;
	quadrant* q;
	WINDOWPLACEMENT wp;

	static HPEN hpen;					// for drawing the zoom box
	static int prev_mouse_x = -1, prev_mouse_y;
	static int dragging = 0, have_box = 0;
	static int allow_mode_change = 1;   // flag indicating whether zoom/pan 
										// mode change is allowed
	static int zoom_mode_pending = 0;   // flag indicating need a change back 
										// to zoom mode
	static int prev_nav_mode = MODE_RTZOOM;
	static int prev_sizing = 0;
	static int prev_max_restore = 0;	// used to detect when window 
										// transitioned from maximized to 
										// restored

	man_calc_struct* m = &main_man_calc_struct;

	switch (nMsg)
	{
	case WM_CREATE: // The window is being created

		hscreen_dc = CreateCompatibleDC(NULL); // screen device context
		hwnd_dialog = CreateDialog(hinstance, MAKEINTRESOURCE(IDD_MAN_DIALOG), 
			hwnd, man_dialog_proc);
		
		// pen for the zoom rectangle; easier to see than PS_DOT
		hpen = CreatePen(PS_SOLID, 2, RGB(0, 0, 0));
		
		setup_sliders();
		set_alg_warning(); // warn if normalized rendering and fast alg

		return FALSE;

	case WM_PAINT: // The window needs to be painted (redrawn).

		hdc = BeginPaint(hwnd, &ps);

		// Blit the mandelbrot bitmap. Could take rectangular regions from 1 to 
		// 4 quadrants.
		// Blits thin horizontal stripes. See comments at striped_blit().
		// Do upper quadrants (UL and UR), then lower (LL and LR).

		// Also optimize for the case where we've recalculated the whole image 
		// (only UL is valid). In this case do a normal blit.
		q = &quad[UL];
		if (q->blit_xsize == m->xsize && q->blit_ysize == m->ysize)
		{
			SelectObject(hscreen_dc, q->handle);
			BitBlt(hdc, 0, 0, q->blit_xsize, q->blit_ysize, hscreen_dc, 0, 0, 
				SRCCOPY);
		}
		else
		{
			striped_blit(&quad[UL], &quad[UR], hdc, hscreen_dc);
			striped_blit(&quad[LL], &quad[LR], hdc, hscreen_dc);
		}
		EndPaint(hwnd, &ps);
		return FALSE;

	case WM_LBUTTONDOWN:

		// Save the coordinates of the mouse cursor.
		mouse_x[0] = LOWORD(lParam);
		mouse_y[0] = HIWORD(lParam);
		mouse_x[1] = LOWORD(lParam); // Init this for switch from zoom to pan also
		mouse_y[1] = HIWORD(lParam);

		// update re/im from any pan offsets and reset offsets
		update_re_im(m, m->pan_xoffs, m->pan_yoffs);
		
		// Get re/im coords at the mouse position, for realtime zoom
		get_mouse_re_im(mouse_x[0], mouse_y[0]);
		
		prev_mouse_x = -1;
		dragging = 1; // dragging either rectangle (for zoom) or image (for pan)

		confine_mouse_cursor();

		if (nav_mode == MODE_PAN)
			SetCursor(hclosed_cursor);	// set closed hand
		if (nav_mode == MODE_RTZOOM)	// if in realtime zoom mode,
			do_rtzoom = RTZOOM_IN;		// set global flag for do_zooming()
		else
			allow_mode_change = 0;		// don't allow mode change while button 
										// is down if not in rtzoom mode

		return FALSE;

	case WM_LBUTTONUP: // Zoom in

		mouse_x[1] = LOWORD(lParam);
		mouse_y[1] = HIWORD(lParam);
		allow_mode_change = 1; // Allow mode change again after button released

		if (dragging)
		{
			dragging = 0;
			if (nav_mode == MODE_ZOOM)
			{
				// Update mandelbrot parms from rectangle and recalculate
				update_re_im_mag(have_box, 1, 
					mouse_x[0], 
					mouse_y[0], 
					mouse_x[1], 
					mouse_y[1]);
				do_man_calculate(1);
			}
			else	// Update image info (excluding iters/sec) after pan or 
					// realtime zoom
				SetWindowText(hwnd_info, get_image_info(0));
		}
		have_box = 0;

		if (zoom_mode_pending) // go back to zoom mode if change was pending
		{
			nav_mode = prev_nav_mode;
			
			// Only if space released in client area (kluge)
			if (GetCursor() != arrow_cursor)
				SetCursor(mag_zoom_cursor);
			zoom_mode_pending = 0;
		}
		// allow dragging during zooming. Possibly some bug here, but I think 
		// I fixed it
		if (nav_mode != MODE_PAN)			// this implements zoom lock
			do_rtzoom = prev_do_rtzoom = 0;	// clear realtime zoom flag for 
											// do_zooming()

		ClipCursor(NULL); // release mouse cursor and capture
		ReleaseCapture();

		return FALSE;

	case WM_MOUSEMOVE:

		mouse_x[1] = LOWORD(lParam);
		mouse_y[1] = HIWORD(lParam);
		
		// Get re/im coords at the mouse position, for realtime zoom
		get_mouse_re_im(mouse_x[1], mouse_y[1]);

		// If user is dragging, draw the zoom rectangle (zoom mode) or pan the 
		// image (pan mode)
		if (nav_mode == MODE_PAN) // panning mode- move the image
		{
			// allow panning using right button drag also
			if (wParam & (MK_LBUTTON | MK_RBUTTON))
			{
				int offs_x, offs_y;

				// Get difference from previous mouse location; use as pan 
				// offset
				offs_x = mouse_x[1] - mouse_x[0];
				offs_y = mouse_y[1] - mouse_y[0];

				mouse_x[0] = mouse_x[1];	// update previous mouse location
				mouse_y[0] = mouse_y[1];

				pan_image(offs_x, offs_y);  // do the pan
			}
		}
		else if (nav_mode == MODE_ZOOM)
		{
			if ((wParam & MK_LBUTTON) && dragging) // zoom rectangle
			{
				hdc = GetDC(hwnd);
				SelectObject(hdc, hpen);

				// not ideal- can be tough to see at times
				SetROP2(hdc, R2_NOTXORPEN);

				// erase previous rectangle, if it exists
				if (prev_mouse_x >= 0 && prev_mouse_x != mouse_x[0])
				{
					Rectangle(hdc, mouse_x[0], mouse_y[0], 
						prev_mouse_x, prev_mouse_y);
					have_box = 1;
				}

				// draw new rectangle
				Rectangle(hdc, mouse_x[0], mouse_y[0],
					prev_mouse_x = mouse_x[1], prev_mouse_y = mouse_y[1]);
				ReleaseDC(hwnd, hdc);
			}
		}
		return FALSE;

	case WM_RBUTTONDOWN: // Zoom out

		mouse_x[0] = LOWORD(lParam);
		mouse_y[0] = HIWORD(lParam);
		mouse_x[1] = LOWORD(lParam); // Init this for switch from zoom to pan also
		mouse_y[1] = HIWORD(lParam);

		// update re/im from any pan offsets and reset offsets
		update_re_im(m, m->pan_xoffs, m->pan_yoffs);
		
		// get re/im coords at the mouse position, for realtime zoom
		get_mouse_re_im(mouse_x[0], mouse_y[0]);

		if (nav_mode == MODE_RTZOOM)	// if in realtime zoom mode,
			do_rtzoom = RTZOOM_OUT;		// set global flag for do_zooming()
		else
			allow_mode_change = 0;		// else don't allow mode change while 
										// button is down

		if (nav_mode == MODE_PAN)
			SetCursor(hclosed_cursor);	// set closed hand

		confine_mouse_cursor();	// also need to confine here, for realtime zoom
		return FALSE;

	case WM_RBUTTONUP: // Zoom out

		mouse_x[1] = LOWORD(lParam);
		mouse_y[1] = HIWORD(lParam);
		allow_mode_change = 1; // allow mode change again after button released

		// Zoom out from current point, and recenter
		if (nav_mode == MODE_ZOOM)
		{
			update_re_im_mag(0, 0, 
				mouse_x[0], 
				mouse_y[0], 
				mouse_x[1], 
				mouse_y[1]);

			do_man_calculate(1);
		}
		if (nav_mode != MODE_PAN)			// this implements zoom lock
			do_rtzoom = prev_do_rtzoom = 0;	// clear realtime zoom flag for 
											// do_zooming()

		ClipCursor(NULL);	// release mouse cursor and capture
		ReleaseCapture();
		return FALSE;

	case WM_MOUSEWHEEL: // Use mousewheel to adjust iterations. (Maybe palette 
						// too, if button down?)

		if (GET_WHEEL_DELTA_WPARAM(wParam) > 0)
			update_iters(1, 0);
		else
			update_iters(0, 1);
		SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
		return FALSE;

	case WM_KEYDOWN:	// Go to pan mode while space is held down (if allowed)

#define PREV_KEYDOWN (1 << 30)

		// aargh... ignore key autorepeats. Were wiping out prev_do_rtzoom.
		if (lParam & PREV_KEYDOWN)
			return TRUE;

		if (allow_mode_change)
		{
			if (wParam == 'Z') // toggle zoom mode
			{
				mag_zoom_cursor = (mag_zoom_cursor == mag_cursor) 
					? rtzoom_cursor : mag_cursor;
				nav_mode = prev_nav_mode = (mag_zoom_cursor == mag_cursor) 
					? MODE_ZOOM : MODE_RTZOOM;

				// Change cursor only if released in client area (kluge)
				if (GetCursor() != arrow_cursor)
					SetCursor(mag_zoom_cursor);
			}
			nav_mode = (wParam == VK_SPACE) ? MODE_PAN : prev_nav_mode;
			if (nav_mode == MODE_PAN)
			{
				mouse_x[0] = mouse_x[1];	// Reset mouse position for pan
				mouse_y[0] = mouse_y[1];
				prev_do_rtzoom = do_rtzoom; // save realtime zoom state
				do_rtzoom = 0;				// stop any realtime zoom
				
				// Fix toggling between hand and arrow in non-client area 
				// (kluge)
				if (GetCursor() == mag_zoom_cursor)
					SetCursor(prev_do_rtzoom ? hclosed_cursor : hopen_cursor);
			}
		}

		// Handle the various hotkeys
		switch (wParam)
		{
		case 'C':	// 'C' toggles the control dialog on and off. Allow to work 
					// in non-fullscreen mode too.

			update_dialog((status ^= STAT_DIALOG_HIDDEN) & 
				STAT_DIALOG_HIDDEN, 0); // 0 = don't move
			break;
		case VK_ESCAPE:	// ESC exits out of fullscreen mode but does not enter 
						// it.
			// deliberate fallthrough if fullscreen
			if (!(status & STAT_FULLSCREEN))
				break;
		case 'F': // 'F' both exits and enters fullscreen mode.
			SendMessage(hwnd_dialog, WM_COMMAND, ID_FULLSCREEN, 0);
			break;
		case 'N': // 'N', 'P', and 'H' do log next/previous and home buttons.
			SendMessage(hwnd_dialog, WM_COMMAND, ID_LOG_NEXT, 0);
			break;
		case 'P':
			SendMessage(hwnd_dialog, WM_COMMAND, ID_LOG_PREV, 0);
			break;
		case 'H':
			SendMessage(hwnd_dialog, WM_COMMAND, ID_HOME, 0);
			break;
		case 'L':	// 'L' toggles palette lock. If locked, palettes in 
					// logfiles are ignored.

			status ^= STAT_PALETTE_LOCKED;
			print_palette_status();
			break;
		case 'I':	// 'I' toggles palette inversion.
			m->pal_xor ^= 0xFFFFFF;
			SendMessage(hwnd_dialog, WM_COMMAND, 
				MAKELONG(IDC_PALETTE, CBN_SELCHANGE), 0);
			print_palette_status();
			break;
		}
		return TRUE;

	case WM_HELP:	// F1 shows help
		SendMessage(hwnd_dialog, WM_COMMAND, ID_HELP_BUTTON, 0);
		return TRUE;

	case WM_KEYUP: // Go back to zoom mode if space released, if allowed.

		if (nav_mode == MODE_PAN)
			if (allow_mode_change)
			{
				// Change cursor only if space released in client area (kluge)
				if (GetCursor() != arrow_cursor)
					SetCursor(mag_zoom_cursor);
				nav_mode = prev_nav_mode;	// Restore old nav mode and
				do_rtzoom = prev_do_rtzoom;	// zooming status
			}
			else
				zoom_mode_pending = 1;	// Else pending: do as soon as mouse 
										// released
		return TRUE;

	case WM_SETCURSOR:	// We get this message whenever the cursor moves in 
						// this window.

		// Set cursor (hand or zoom) based on nav mode. Also set keyboard focus 
		// to this window for key detection. Eliminates need to click window 
		// first to set focus
		SetFocus(hwnd);

		// only set zoom/hand cursor in client area
		if (LOWORD(lParam) == HTCLIENT)
		{
			SetCursor(nav_mode == MODE_PAN ? hopen_cursor : mag_zoom_cursor);
			return TRUE;
		}
		break; // let system set cursor outside client area (resize, arrow, etc)

		// This message comes after the user finishes a drag or movement, but, 
		// annoyingly, not after a maximize or restore. There's some extra code 
		// to handle those in WM_WINDOWPOSCHANGED

	case WM_EXITSIZEMOVE:
		// if previous resize/move was a resize, 
		if (prev_sizing)
			// recalculate image if enabled
			if (cfg_default.options.c_val & OPT_RECALC_ON_RESIZE)
			{
				// If we recalculate here, a partial bad-state window appears 
				// while the calculation is going on. Recalculate in main loop 
				// instead.
				status |= STAT_RECALC_IMMEDIATELY;
			}
		return FALSE;

		// Called on window sizing or changing position. Pretty wasteful to 
		// call create_bitmap (which frees and reallocates all arrays with 
		// every pixel of movement) here. Memory fragmentation? Changing this 
		// currently causes some issues. Maybe fix later.

		// Seems like we always get this at startup before the paint message.

	case WM_WINDOWPOSCHANGED:
		GetWindowPlacement(hwnd, &wp);
		// only do this if not minimizing window
		if (wp.showCmd != SW_SHOWMINIMIZED)
		{
			// Move dialog box along with main window
			update_dialog(status & STAT_DIALOG_HIDDEN, 1); // 1 = move

			GetClientRect(hwnd, &rc);	// calculate new mandelbrot image size
			m->xsize = rc.right - rc.left;
			m->ysize = rc.bottom - rc.top;

			if (m->xsize < MIN_SIZE)	// clip min size
				m->xsize = MIN_SIZE;
			if (m->ysize < MIN_SIZE)
				m->ysize = MIN_SIZE;

			// Create arrays and bitmaps used in calculations (will 
			// deallocate/resize as needed- returns > 0 if resized). Doesn't do 
			// anything if size didn't change. Set prev_sizing so 
			// WM_EXITSIZEMOVE can know whether previous op was a resize.

			// rename this- now does a lot more than create a bitmap
			if (prev_sizing = create_bitmap(m->xsize, m->ysize))
			{
				// update size info if resized
				SetWindowText(hwnd_info, get_image_info(0));

				// update sizes for image save
				SetDlgItemInt(hwnd_dialog, IDC_SAVE_XSIZE, m->xsize, FALSE);
				SetDlgItemInt(hwnd_dialog, IDC_SAVE_YSIZE, m->ysize, FALSE);
			}

			// Send a WM_EXITSIZEMOVE message if window was maximized/restored
			if (wp.showCmd != prev_max_restore)
			{
				SendMessage(hwnd, WM_EXITSIZEMOVE, 0, 0);
				prev_max_restore = wp.showCmd;
			}
			// Save main window rect for fullscreen mode if not maximized and 
			// not already in fullscreen
			if (wp.showCmd != SW_SHOWMAXIMIZED && !(status & STAT_FULLSCREEN))
				GetWindowRect(hwnd, &main_rect);
			return FALSE;
		}
		return TRUE;

	case WM_COMMAND:
		return FALSE;

	case WM_DESTROY:	// The window is being destroyed, close the application
		PostQuitMessage(0);
		return FALSE;
	}
	// If we don't handle a message completely we hand it to the 
	// system-provided default window function.
	return DefWindowProc(hwnd, nMsg, wParam, lParam);
}

void fancy_intro(void)
{
#define MAG_STEP 1.07 // slow down a bit from previous versions

	man_calc_struct* m = &main_man_calc_struct;

	set_home_image();
	m->max_iters = 64;
	m->mag = MAG_START;
	do_rtzoom = 1; // prevent status line from being updated
	do
	{
		do_man_calculate(1);
		m->mag *= MAG_STEP;
	} while (m->mag <= 1.35);

	set_home_image();
	do_man_calculate(1);
	do_rtzoom = 0;
	// "resized" initially, but no need to calc again
	status &= ~STAT_RECALC_IMMEDIATELY;

	SetWindowText(hwnd_info, get_image_info(1)); // print first info, status
	print_status_line(0);

	file_tot_time = 0.0; // don't count intro time in file total time
}

int CALLBACK 
WinMain(HINSTANCE hInst, HINSTANCE hPrev, LPSTR lpCmd, int nShow)
{
	MSG msg;
	WNDCLASSEX wndclass;
	static char* classname = "ManWin";
	man_calc_struct* m = &main_man_calc_struct;
	man_calc_struct* s = &save_man_calc_struct;

	hinstance = hInst;

	// Read any default settings from quickman.cfg. Do this early, because it 
	// can change almost anything used below.
	read_cfg_file();
	get_cpu_info();
	get_system_metrics();
	man_init();
	if (!(num_builtin_palettes = init_palettes(DIVERGED_THRESH, m, s)))
		return 0;

	// create a window class for our main window
	memset(&wndclass, 0, sizeof(WNDCLASSEX));
	wndclass.lpszClassName = classname;
	wndclass.cbSize = sizeof(WNDCLASSEX);
	wndclass.style = CS_HREDRAW | CS_VREDRAW;
	wndclass.lpfnWndProc = MainWndProc;
	wndclass.hInstance = hInst;
	// use mini mandelbrot icon
	wndclass.hIcon = LoadIcon(hInst, MAKEINTRESOURCE(IDI_MAN));
	
	wndclass.hIconSm = NULL;	// derive small icon from normal one
	// we set other cursors based on nav mode
	wndclass.hCursor = arrow_cursor = LoadCursor(NULL, IDC_ARROW);
	
	wndclass.hbrBackground = NULL;	// don't need this
	RegisterClassEx(&wndclass);

	// Load cursors - wait, zoom, hand open, and hand closed
	wait_cursor = LoadCursor(NULL, IDC_WAIT);
	mag_cursor = LoadCursor(hInst, MAKEINTRESOURCE(IDC_MAG));
	rtzoom_cursor = mag_zoom_cursor = LoadCursor(hInst, 
		MAKEINTRESOURCE(IDC_RTZOOM));
	hopen_cursor = LoadCursor(hInst, MAKEINTRESOURCE(IDC_HAND_OPEN));
	hclosed_cursor = LoadCursor(hInst, MAKEINTRESOURCE(IDC_HAND_CLOSED));

	// Create our main window
	hwnd_main = CreateWindow(
		classname,
		"QuickMAN 1.10  |  F1: Help",
		WS_OVERLAPPEDWINDOW,	// Style
		140, // CW_USEDEFAULT,	// Initial x; come up at a good location on my 
								// laptop...
		20,  // CW_USEDEFAULT,	// Initial y
		m->xsize + x_border,	// Size
		m->ysize + y_border,
		NULL,					// No parent window
		NULL,					// No menu
		hInst,					// This program instance
		NULL					// Creation parameters
	);

	// Init common controls; required for < WinXP, and apparently for some 
	// Pentium 4 systems
	InitCommonControls();

#ifndef USE_PERFORMANCE_COUNTER
	// Set timer resolution to 1ms (if supported). Seems to take awhile to take 
	// effect. timegettime() will return inaccurate results if it's called too 
	// soon after this...
	timeBeginPeriod(1); // fancy_intro will absorb initialization time
#endif

	UpdateWindow(hwnd_main);
	update_dialog(1, 1);          // 1 = hide, 1 = move
	ShowWindow(hwnd_main, nShow); // causes dialog to become visible
	UpdateWindow(hwnd_dialog);

	// Add user palettes and logfiles to their dropdown menus
	add_user_palettes_and_logfiles();
	
	// Read logfile - add_user_palettes() must be called before this
	log_read(logfile, 
		"\nDid you extract all the files from the QuickMAN .zip archive?", 1);
	fancy_intro(); // Zoom in to home image

	while (1) // Main loop
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
				break;
			
			// only process messages not for dialog box
			if (!IsDialogMessage(hwnd_dialog, &msg))
			{
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}
		else
		{
			// Do heavy computation here
			if (!do_zooming() && !do_panning() && !do_recalc())
				Sleep(2);	// don't use 100% of CPU when idle. Also see 
							// do_panning()
		}
	}

#ifndef USE_PERFORMANCE_COUNTER
	timeEndPeriod(1);
#endif

	free_man_mem(&main_man_calc_struct);
	free_man_mem(&save_man_calc_struct);

	return (int)msg.wParam;
}