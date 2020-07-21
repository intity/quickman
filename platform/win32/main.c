/******************************************************************************
 * main.c
 * -> Win32 GUI for the QuickMAN SSE/SSE2-based Mandelbrot Set calculator
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
#define ALIGNED_(x) __declspec(align(x))

#include <windows.h>
#include <mmsystem.h>
#include <process.h>
#include <commctrl.h>

#ifdef _DEBUG
#include <crtdbg.h>
#endif

#include "stdafx.h"
#include "resource.h"
#include "..\..\src\quickman.h"

void do_man_calculate(int recalc_all);

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

/**
 * Combo box initialization strings/defines. String order should correspond
 * to the above define order for PRECISION_* and ALG_*.
 */
static char* precision_strs[] = { "Auto", "Single", "Double", "Extended" };

/**
 * Combo box initialization: algorithm list.
 */
static char* algorithm_strs[] = {
	"Fast, AMD", 
	"Exact, AMD", 
	"Fast, Intel", 
	"Exact, Intel", 
	"Fast, C", 
	"Exact, C"
};

/**
 * Combo box initialization: 'Standard' or 'Normalized' iteration count.
 */
static char* rendering_strs[] = { "Standard", "Normalized" };

/**
 * Combo box initialization: not all of these will be used.
 */
static char* threads_strs[] = {
	"1", "2", "4", "8", "16", "32", "64", "128", "256"
};

/* ---------------------- Global mandelbrot parameters --------------------- */

static int prev_xsize;				// previous sizes, for restoring window
static int prev_ysize;
static double mouse_re;				// re/im coordinates of the mouse position
static double mouse_im;
static double zoom_start_mag;		// starting magnification, for zoom button
static unsigned int num_palettes;	// total number of builtin palettes

static char log_file[256] = LOG_FILE;	// default logfile
static char img_file[256] = IMG_FILE;	// default save filename
static char pal_file[256];			// filename of current user palette file

static int nav_mode = NAV_RTZOOM;	// navigation mode
static int img_save = 0;	// used in from do_save for call man_calculate
static int do_rtzoom = 0;	// if nonzero, do a realtime zoom in/out
static int prev_do_rtzoom;	// previous state of do_rtzoom

/**
 * Presets for the logfile combo box. Only need to add files here that we want
 * to be at the top of the box. Others will be added from the current directory.
 */
char* file_strs[] = { log_file, "auto_panzoom.log" };

/**
 * For zoom button benchmarking (helps measure overhead).
 */
static double zoom_start_time;

/**
 * Mouse position: index 0 is position on initial button press; index 1 is 
 * current position.
 */
static int mouse_x[2], mouse_y[2];

// GUI stuff
static HCURSOR rtzoom_cursor;
static HCURSOR hopen_cursor;
static HCURSOR hclosed_cursor;
static HCURSOR arrow_cursor;
static HCURSOR wait_cursor;
static HCURSOR mag_cursor;
static HCURSOR mag_zoom_cursor;		// either mag or rtzoom, depending on mode
static RECT main_rect;				// rectangle for the main window
static HWND hwnd_main;				// handle for the main window
static HWND hwnd_dialog = 0;		// handle for the dialog box
static HWND hwnd_info;				// handle for the info text areas
static HWND hwnd_status_1;			// handle for the status text areas
static HWND hwnd_status_2;
static HWND hwnd_iters;				// and other controls
static HINSTANCE hinstance = 0;		// the application
static HDC hscreen_dc = NULL;		// screen device context

// variables for calculating window sizes and such
static int x_border;
static int y_border;

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

// used for normal calculation
ALIGNED_(64) man_calc_struct main_man_calc_struct;

// used for saving images
ALIGNED_(64) man_calc_struct save_man_calc_struct;

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

#define PAN_KEY (KEY_RIGHT | KEY_LEFT | KEY_UP | KEY_DOWN) // any pan key

/**
 * Some magic constants
 */
#define PAN_STEP_DIV		150000.0
#define OVERHEAD_FACTOR		100000

/**
 * These adjust the pan filter constant based on image size (filter constant
 * now comes from config file- adds acceleration and deceleration to movements).
 */
#define PFC_SLOPE_FACTOR   (1600.0 * 1140.0 - 700.0 * 700.0)
#define PFC_OFFS_FACTOR    (700.0 * 700.0)

static double cur_pan_xstep = 0.0;
static double cur_pan_ystep = 0.0;
static double pan_xstep_accum = 0.0;
static double pan_ystep_accum = 0.0;

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
static settings cfg = {
	/* name				c_val			d_val			min		max */
	{ "panrate",		DEF_PAN_RATE,	DEF_PAN_RATE,	0,		MAX_PAN_RATE	},
	{ "pan",			0,				0,				0,		0xFFFF			},
	{ "zoomrate",		DEF_ZOOM_RATE,	DEF_ZOOM_RATE,	0,		MAX_ZOOM_RATE	},
	{ "zoom",			0,				0,				0,		0xFFFF			},
	{ "xsize",			700,			700,			0,		0xFFFF			},
	{ "ysize",			700,			700,			0,		0xFFFF			},
	{ "maxiters_color", 0,				0,				0,		0xFFFFFF		},
	{ "pal_xor",		0,				0,				0,		0xFFFFFF		},
	{ "options",		OPT_DEFAULT,	OPT_DEFAULT,	0,		0xFFFF			},
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
static settings cfg_cur;

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
	for (i = 0; i < CFG_SIZE; i++)
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
	for (i = 0; i < CFG_SIZE; i++)
		if (SETTING_AUTORESET(&d[i]))
			d[i].c_val = d[i].d_val;
}

/**
 * Set all val fields to -1, which invalidates all settings.
 */
void invalidate_settings(settings* dest)
{
	int i;
	setting* d;
	d = (setting*)dest; // treat struct as array
	for (i = 0; i < CFG_SIZE; i++)
		d[i].c_val = -1;
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

	if (fopen_s(&fp, file, bin ? "rb" : "rt"))
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
	man_calc_struct* m = &main_man_calc_struct;

	// Initialize cur file settings structure to all invalid (no change)
	invalidate_settings(&cfg_cur);

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

		s = (setting*)&cfg;		// treat structs as arrays
		f = (setting*)&cfg_cur;
		for (j = 0; j < sizeof(cfg) / sizeof(setting); j++)
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

		entry->log_settings = cfg_cur; // Copy any settings found above

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
				num_palettes - 1,
				(LPARAM)&str[j + 1]);

			entry->palette = (i != CB_ERR) ? i : DEFAULT_PAL;
		}
	}
	return 1;
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
	man_calc_struct* m = &main_man_calc_struct;

	log_count = 0;
	if (init_pos)
	{
		log_pos = -1;
		m->calc_time = 0.0; // for benchmarking
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
	if (m->palette < num_palettes)
		sprintf_s(p, sizeof(p), "%d", m->palette);
	else
		sprintf_s(p, sizeof(p), "\"%s\"", pal_file);

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
	if (!(m->status & STAT_PALETTE_LOCKED))
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

	invalidate_settings(&cfg_cur); // initialize all to "no change"

	// NULL = no error message on failure to read
	log_read(CFG_FILE, NULL, 1);

	// 1 = copy to default_val also
	copy_changed_settings(&cfg, &cfg_cur, 1);

	// maybe eliminate these separate variables later
	m->xsize = prev_xsize = cfg.xsize.c_val;
	m->ysize = prev_ysize = cfg.ysize.c_val;

	// just so 1st dialog box doesn't get div by 0
	m->min_dimension = m->xsize > m->ysize ? m->ysize : m->xsize;
	m->max_iters_color = cfg.max_iters_color.c_val;
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
 * Print frames/sec info and current algorithm indicator. Current algorithm can
 * temporarily change from Fast to Exact while panning, to improve performance.
 * See man_calculate_threaded.
 */
void print_fps_status_line(double fps, double avg_fps, double eff)
{
	char s[256];
	man_calc_struct* m = &main_man_calc_struct;

	// if currently saving, keep saving status line
	if (m->status & STAT_DOING_SAVE)
		return;

	// Interval and average frames/sec. Removed "AVG" so large frame rates 
	// don't get cut off
	sprintf_s(s, sizeof(s), "%c Fps %3.0f/%-3.0f", m->cur_alg & ALG_EXACT 
		? 'E' : 'F', fps, avg_fps);
	SetWindowText(hwnd_status_1, s);

	// Iteration percentage: mandelbrot calculation time / total time
	sprintf_s(s, sizeof(s), "Iter %2.0f%%", eff);
	SetWindowText(hwnd_status_2, s);
}

/**
 * Print the text in the INFO area of the toolbar.
 */
void print_image_info(int update_iters_sec)
{
	int i;
	char s[1024 + 32 * MAX_THREADS], tmp[256];
	man_calc_struct* m = &main_man_calc_struct;
	thread_state* t;

	get_image_info(m, update_iters_sec);

	sprintf_s(s, sizeof(s), 
		"Real\t%-16.16lf\r\n"
		"Imag\t % -16.16lf\r\n"
		"Mag\t%-16lf\r\n"
		"\r\n"
		"Size\t%u x %u\r\n"
		"Time\t%-3.3fs\r\n"
		"Iters/s\t%-4.4gM (%-.2f GFlops)\r\n"
		"\r\n"
		"Avg iters/pixel\t%-.1lf\r\n"
		"Points guessed\t%-.1lf%%\r\n"
		"Total iterations\t%-.0lf\r\n",

		m->re + get_re_im_offs(m, m->pan_xoffs),
		m->im - get_re_im_offs(m, m->pan_yoffs),
		m->mag,
		m->xsize,
		m->ysize,
		m->iter_time,
		m->ictr_m,
		m->ictr_g,
		m->avg_iters,
		m->guessed_pct,
		(double)m->ictr
	);

	sprintf_s(tmp, sizeof(tmp), "\r\nThread load %%\tCur\tTotal\r\n");
	strcat_s(s, sizeof(s), tmp);

	for (i = 0; i < m->num_threads; i++)
	{
		t = &m->thread_states[i];

		sprintf_s(tmp, sizeof(tmp), 
			"Thread %d\t%#3.3g\t%#3.3g\r\n", i, 
			t->cur_pct,
			t->tot_pct);
		strcat_s(s, sizeof(s), tmp);
	}

	sprintf_s(tmp, sizeof(tmp),
		"Efficiency %%\t%#3.3g\t%#3.3g\r\n"
		"\r\n"
		"Total calculate time:\t%-.3lfs",
		m->max_cur_pct,
		m->max_tot_pct,
		m->calc_time);

	strcat_s(s, sizeof(s), tmp);

	SetWindowText(hwnd_info, s);
}

/**
 * Get frames/sec during one interval, average frames per sec, and iteration
 * time percentage (mandelbrot calculation time / pan or zoom time). Assumes
 * iter_time is set to the mandelbrot calculation time. Op_time should be the
 * total time for the pan or zoom operation. Call once per frame.
 */
void update_benchmarks(double op_time, int update_iters_sec)
{
#define UPDATE_INTERVAL_TIME	0.25 // update about 4 times per sec

	double fps, avg_fps, eff;
	man_calc_struct* m = &main_man_calc_struct;

	m->interval_frames++;
	m->total_frames++;
	m->interval_time += op_time; // Update interval and total operation times
	m->total_time += op_time;

	// Update time spent only on mandelbrot calculation
	m->calc_total_time += m->iter_time;

	// Update time spent only on mandelbrot calculation this interval
	m->calc_interval_time += m->iter_time;

	// Update status line every interval
	if (m->interval_time >= UPDATE_INTERVAL_TIME)
	{
		fps = (double)m->interval_frames / m->interval_time;
		avg_fps = (double)m->total_frames / m->total_time;

		// interval iteration %
		eff = 100.0 * m->calc_interval_time / m->interval_time;

		// average iteration %
		// eff = 100.0 * calc_total_time / total_time;

		print_fps_status_line(fps, avg_fps, eff);
		print_image_info(update_iters_sec);	// update image info

		m->interval_frames    = 0;
		m->interval_time      = 0.0;
		m->calc_interval_time = 0.0;
	}
}

/**
 * Reset the pan filter state and accumulators.
 */
void reset_pan_state(man_calc_struct* m)
{
	cur_pan_xstep	= 0.0;
	cur_pan_ystep	= 0.0;
	pan_xstep_accum = 0.0;
	pan_ystep_accum = 0.0;
}

/**
 * This calculates the x and y steps for panning based on keys pressed.
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
	// not in pan lock mode: see if we need to go into it
	else if ((pkey & KEY_SHIFT) && (pkey & PAN_KEY))
	{
		key_lock = pkey & (PAN_KEY | KEY_CTRL);
		wait_release = 1;
	}

	// Get the step scale from the pan rate slider
	pan_step_scale = pan_step_scales[cfg.pan_rate.c_val];

	// Adjust the step based on the image size, to try to compensate for 
	// size-based variable frame rates
	pan_step = pan_step_scale * (double)
		(m->img_size + OVERHEAD_FACTOR) * (1.0 / PAN_STEP_DIV);

	// Pan filter settings now come from config file. The 0.0001 scales large 
	// config file integers to the values needed here.
	pfcmin = 0.0001 * (double)cfg.pfcmin.c_val;
	pfcmax = 0.0001 * (double)cfg.pfcmax.c_val;
	pfc_slope = (pfcmax - pfcmin) * (1.0 / PFC_SLOPE_FACTOR);
	pfc_offs = pfcmin - pfc_slope * PFC_OFFS_FACTOR;

	pan_filter_const = (double)m->img_size * pfc_slope + pfc_offs;
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
	cur_pan_xstep *= (tmp = 1.0 - pan_filter_const);
	cur_pan_xstep += xs * pan_filter_const;
	cur_pan_ystep *= tmp;
	cur_pan_ystep += ys * pan_filter_const;

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

	// else decrement and pulse on transition to 0
	else if (stopped_counter && !--stopped_counter)
	{
		stopped = 1;

		// Clear out old data. Without this, there can be little artifacts 
		// when restarting.
		reset_pan_state(m);
		pulse = 1;
	}

	return pulse;
}

/**
 * Get re/im coordinates at the mouse position (mx, my), for realtime zoom.
 */
void get_mouse_re_im(int mx, int my)
{
	man_calc_struct* m = &main_man_calc_struct;

	mx -= (m->xsize >> 1); // Get offset from image center
	my -= (m->ysize >> 1);

	mouse_re = m->re + get_re_im_offs(m, mx);
	mouse_im = m->im - get_re_im_offs(m, my);
}

/**
 * Pan the image using the keyboard. Super cool...
 *
 * Returns 1 if it did a pan, else 0 (if idle).
 */
int do_panning()
{
	int xstep = 0, ystep = 0;
	static double start_time;
	static double pan_time = -1.0;
	man_calc_struct* m = &main_man_calc_struct;

	// Update coords, xstep, and ystep based on keys pressed. Returns a 1-clock 
	// pulse when panning transitioned from active to stopped.

	if (get_pan_steps(&xstep, &ystep, 0))
	{
		print_image_info(0);	// update image info

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

		pan_image(m, xstep, ystep);
		do_man_calculate(0);

		// get time since last screen update
		pan_time = get_seconds_elapsed(start_time);

		// skip update if the whole image was recalculated (due to size change, 
		// etc). messes up the averages.
		if (!m->all_recalculated)
		{
			// update the fps, average fps, and iteration %
			update_benchmarks(pan_time, 0);
		}
		start_time = get_timer();	// restart timer for next update

		return 1;
	}
	return 0;
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
int do_zooming()
{
	double start_time;
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
			if (!(do_rtzoom & RTZOOM_BTN))
				prev_do_rtzoom = do_rtzoom;
			else
				prev_do_rtzoom = 0;

		do_rtzoom = 0;
	}
	else if (prev_do_rtzoom)
	{
		// update re/im from any pan offsets and reset offsets
		update_re_im(m, m->pan_xoffs, m->pan_yoffs);

		// update mouse coords after any pan
		get_mouse_re_im(mouse_x[1], mouse_y[1]);

		// reset pan so we don't get extra movement
		reset_pan_state(m);

		// after stopping zoom
		do_rtzoom = prev_do_rtzoom;
	}

	if (!do_rtzoom) // Return if not zooming
		return 0;

	// update re/im from any pan offsets and reset offsets
	update_re_im(m, m->pan_xoffs, m->pan_yoffs);

	step = rtzoom_mag_steps[cfg.zoom_rate.c_val];

	start_time = get_timer();

	if (do_rtzoom & RTZOOM_IN)
		m->mag *= step;
	else
	{
		m->mag /= step;
		if (m->mag < MAG_MIN)
			m->mag = MAG_MIN;
	}

	if (!(do_rtzoom & RTZOOM_BTN)) // if zooming using the mouse
	{
		// Set the new image center re/im to keep the position at mouse[1]
		// at the same point on the screen.

		mx = mouse_x[1] - (m->xsize >> 1); // Get offset from image center
		my = mouse_y[1] - (m->ysize >> 1);

		m->re = mouse_re - get_re_im_offs(m, mx);
		m->im = mouse_im + get_re_im_offs(m, my);
	}
	// if zooming using the button, stop when we hit the start mag
	else if (m->mag > zoom_start_mag)
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
		m->calc_time = get_seconds_elapsed(zoom_start_time);

		print_image_info(1);	// update image info
	}
	return 1;
}

/**
 * Simple function for recalculating after a window size change (if enabled).
 */
int do_recalc()
{
	man_calc_struct* m = &main_man_calc_struct;

	if (m->status & STAT_RECALC_IMMEDIATELY)
	{
		do_man_calculate(1);
		m->status &= ~STAT_RECALC_IMMEDIATELY;
	}
	return 0;
}

/* --------------------------- Queuing functions --------------------------- */

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
	man_calc_struct* m;
	thread_state* t;
	stripe* s;

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
		xend   = s->xend;
		ystart = s->ystart;
		yend   = s->yend;

		// Optimization for panning: set alg to exact mode for very thin 
		// regions. Due to the Fast algorithm's 4x4 cell size it often computes 
		// more pixels than Exact for these regions. Effect is most apparent 
		// with high iter count images.

#define FE_SWITCHOVER_THRESH  2 // only do it for 1-pixel wide regions for now. 
								// best value TBD...

		m->cur_alg = m->algorithm;
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
			points_guessed = fast_wave_alg(m, ps_ptr, s);
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
		m->queue_point(m, ps_ptr, m->iter_data + m->img_size);

	// Thread 0 always runs in the master thread, so doesn't need to signal. 
	// Save overhead.
	if (t->thread_num)
		SetEvent(t->done_event);	// For other threads, signal master thread 
									// that we're done

	return 0;
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
static double man_calculate(int xstart, int xend, int ystart, int yend)
{
	double start_time;
	int i, xsize, ysize, step, thread_ind, stripe_ind, num_stripes, frac;
	int frac_step, this_step;
	stripe* s = NULL;
	man_calc_struct* m;
	
	m = img_save ? &save_man_calc_struct : &main_man_calc_struct;
	m->all_recalculated = 0;

	// if need recalculation, recalculate all. No effect for saving
	if (m->status & STAT_NEED_RECALC)
	{
		xstart = 0;	// reset rectangle to full screen
		xend = m->xsize - 1;
		ystart = 0;
		yend = m->ysize - 1;
		m->status &= ~STAT_NEED_RECALC;
		m->all_recalculated = 1;
	}

	man_setup(m, xstart, xend, ystart, yend);

	xsize = xend - xstart + 1;
	ysize = yend - ystart + 1;

	num_stripes = m->num_stripes;

	// Need to check min/max here (couldn't be checked automatically by 
	// log-reading function)
	if (num_stripes < 1)
		num_stripes = 1;
	if (num_stripes > MAX_STRIPES)
		num_stripes = MAX_STRIPES;

	// Now multiply by num_threads to get total number of stripes for the 
	// image.
	num_stripes <<= m->threads_i;

	// With pathologically small images, some threads may not calculate 
	// anything.
	for (i = 0; i < m->num_threads; i++)
		m->thread_states[i].num_stripes = 0;

	// Start at the last thread, so that thread 0 gets any leftovers at the end. 
	// Thread 0 is the master and doesn't suffer the overhead of being spawned, 
	// so it should get the extra work.

	thread_ind = m->num_threads - 1;
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
		// Dual Opteron 280, Double, Fast, 4 threads, 4 stripes/thread, 
		// tune.log:
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
				thread_ind = m->num_threads - 1;
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
				thread_ind = m->num_threads - 1;
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

	for (i = 1; i < m->num_threads; i++)
		QueueUserWorkItem(man_calculate_threaded, &m->thread_states[i],
			WT_EXECUTELONGFUNCTION | (MAX_QUEUE_THREADS << 16));

	// Could also queue thread 0 too, then have this master thread display the 
	// progress if the calculation is really slow...
	man_calculate_threaded(&m->thread_states[0]);

	if (m->num_threads > 1)
	{
		WaitForMultipleObjects(m->num_threads - 1,
			&m->thread_done_events[1],
			TRUE, INFINITE); // wait till all threads are done
	}

	if (!(m->flags & FLAG_IS_SAVE)) // don't update these if doing save
	{
		m->calc_time += get_seconds_elapsed(start_time);
		return m->calc_time;
	}
	else
		return 0.0;
}

/* --------------------------- Palette functions --------------------------- */

/**
 * Create a palette from a BMP file. Cool...
 *
 * Only the last horizontal line of the BMP is used. Thus the palette size is
 * the width of the BMP image. Uses the same user palette array that's used for
 * text palettes.
 */
int load_palette_from_bmp(FILE* fp, man_calc_struct* m)
{
	BITMAPFILEHEADER head;
	BITMAPINFOHEADER info;
	unsigned int n;

	// Set user palette to invalid by default (if valid, num_valid_palettes 
	// will be NUM_PALETTES + 1)
	m->num_valid_palettes = NUM_PALETTES;

	// Make sure we read good header and info structures and can seek to the 
	// start of the bitmap data (l->r evaluation order guaranteed). Must have 
	// signature "BM" and be an uncompressed 24-bit bitmap.

	if (!fread(&head, sizeof(head), 1, fp) || head.bfType != 0x4D42 || // "BM"
		!fread(&info, sizeof(info), 1, fp) || info.biSize != sizeof(info) ||
		info.biPlanes != 1 || 
		info.biCompression != BI_RGB ||
		info.biBitCount != 24 ||
		fseek(fp, head.bfOffBits, SEEK_SET))
		return 0;

	// Read in WIDTH 24-bit values. Limit any insane palette sizes from 
	// corrupted files
	if ((n = info.biWidth) > MAX_PALETTE_SIZE)
		return 0;

	return load_palette_from_array(fp, m, n);
}

/**
 * New threaded palette mapping function. Called from apply_palette (the 
 * interface to the rest of the code).
 */
unsigned int CALLBACK apply_palette_threaded(pal_work* p)
{
	man_calc_struct* m = (man_calc_struct*)p->calc_struct;

	calc_palette_for_thread(p);

	// Tell the master thread we're done. Thread 0 is the master thread, so 
	// doesn't need to signal
	if (p->thread_num > 0)
		SetEvent(m->pal_events[p->thread_num]);
	return 0;
}

void apply_palette(unsigned int* dest, unsigned int* src,
	unsigned int xsize,
	unsigned int ysize)
{
	unsigned int i, nt;
	man_calc_struct* m;

	m = img_save ? &save_man_calc_struct : &main_man_calc_struct;

	nt = calc_palette(m, dest, src, xsize, ysize);

	// If more than one thread, spawn new threads
	for (i = 1; i < nt; i++)
		QueueUserWorkItem(apply_palette_threaded, &m->pal_work_array[i],
			WT_EXECUTELONGFUNCTION | (MAX_QUEUE_THREADS << 16));

	// Do some (or all) of the work here in the master thread
	apply_palette_threaded(&m->pal_work_array[0]);

	// Wait till all threads are done
	if (nt > 1)
		WaitForMultipleObjects(nt - 1, &m->pal_events[1], TRUE, INFINITE);
}

/* -------------------------- GUI/misc functions --------------------------- */

/**
 * Get the number of stripes per thread based on number of threads (extract 
 * bits 3-0 for 1 thread, 7-4 for two threads, etc).
 */
int get_num_stripes(int threads_i)
{
	return (cfg.stripes.c_val >> (threads_i << 2)) & 0xF;
}

/**
 * Find a string in a strings array and return its index, or -1 if not found.
 */
int get_string_index(char* str, char** strs, int num_strs)
{
	int i;
	for (i = 0; i < num_strs; i++)
		if (!strcmp(str, strs[i]))
			return i;
	return -1;
}

/**
 * Detect whether the CPU supports SSE2 and conditional move instructions; used
 * to set algorithms. Also detect the number of cores.
 */
void get_cpu_info()
{
	unsigned int vendor[4];
	unsigned int features;
	SYSTEM_INFO info;
	man_calc_struct* m = &main_man_calc_struct;

	vendor[3] = 0;

#define FEATURE_SSE		0x02000000
#define FEATURE_SSE2	0x04000000
#define FEATURE_CMOV	0x00008000

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
		m->algorithm = ALG_FAST_ASM_AMD;
	else
		m->algorithm = ALG_FAST_ASM_INTEL;

	// Set exact alg if configured by quickman.cfg options field
	if (cfg.options.c_val & OPT_EXACT_ALG)
		m->algorithm |= ALG_EXACT;

	m->sse_support = 0;
	if ((features & (FEATURE_SSE | FEATURE_CMOV)) == (FEATURE_SSE | FEATURE_CMOV))
		m->sse_support = 1;
	if ((features & (FEATURE_SSE2 | FEATURE_CMOV)) == (FEATURE_SSE2 | FEATURE_CMOV))
		m->sse_support = 2;

	if (m->sse_support < 2)
	{
		MessageBox(NULL,
			"Your (obsolete) CPU does not support SSE2 instructions.\r\n"
			"Performance will be suboptimal.", "Warning",
			MB_OK | MB_ICONSTOP | MB_TASKMODAL);

		// Ok to stay in auto precision mode with only sse support- will switch
		// to C algorithm for double. If no sse support, no choice but to use C
		if (!m->sse_support)
			m->algorithm = ALG_FAST_C;
	}

	// Set the default number of threads to the number of cores. Does this 
	// count a hyperthreading single core as more than one core? Should ignore 
	// these as hyperthreading won't help.
	GetSystemInfo(&info);
	m->num_threads = info.dwNumberOfProcessors;

	// Convert number of threads (cores) to a selection index for the dropdown 
	// box

	// Get log2(num_threads)
	for (m->threads_i = 0; m->threads_i <= MAX_THREADS_IND; m->threads_i++)
		if ((1 << m->threads_i) >= m->num_threads)
			break;

	m->num_threads = 1 << m->threads_i;
	m->num_stripes = get_num_stripes(m->threads_i);
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
	int i, err;
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
			DeleteObject(m->quad[i].handle);
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
		m->quad[i].handle = CreateDIBSection
		(NULL, &bmi, DIB_RGB_COLORS, (void**)&m->quad[i].bitmap_data, NULL, 0);

		if (!m->quad[i].handle)
			err = 1;
	}

	if (!hscreen_dc || err || !alloc_man_mem(m, width, height))
	{
		MessageBox(NULL, "Error allocating storage arrays.",
			NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
		// really should exit cleanly here... subsequent code will crash
	}

	// update re/im from any pan offsets and reset offsets resized; need to 
	// recalculate entire image
	update_re_im(m, m->pan_xoffs, m->pan_yoffs);
	m->status |= STAT_NEED_RECALC;

	prev_width = width;
	prev_height = height;

	// set smaller dimension
	m->min_dimension = (width < height) ? width : height;

	set_wave_ptr_offs(m);	// set wave pointer offsets

	reset_quadrants(m);		// reset to recalculate all
	reset_fps_values(m);	// reset frames/sec timing values
	reset_pan_state(m);		// reset pan filters and movement state

	return 1;
}

/**
 * Initialize the ComboBox control.
 */
void combobox_init(HWND hwnd, int idc)
{
	int i, seln, size;
	char** strs;
	man_calc_struct* m = &main_man_calc_struct;

	switch (idc)
	{
	case IDC_PALETTE:
		seln = m->palette;
		strs = palette_strs;
		size = NUM_ELEM(palette_strs);
		break;
	case IDC_RENDERING:
		seln = m->rendering;
		strs = rendering_strs;
		size = NUM_ELEM(rendering_strs);
		break;
	case IDC_PRECISION:
		seln = m->precision;
		strs = precision_strs;
		size = NUM_ELEM(precision_strs);
		break;
	case IDC_ALGORITHM:
		seln = m->algorithm;
		strs = algorithm_strs;
		size = NUM_ELEM(algorithm_strs);
		break;
	case IDC_THREADS:
		seln = m->threads_i;
		strs = threads_strs;
		size = MAX_THREADS_IND + 1;
		break;
	case IDC_LOGFILE:
		seln = 0;
		strs = file_strs;
		size = NUM_ELEM(file_strs);
		break;
	}

	// Initialize all the strings
	for (i = 0; i < size; i++)
		SendDlgItemMessage(hwnd, idc, CB_ADDSTRING, 0, (LPARAM)strs[i]);

	// Set the default selection
	SendDlgItemMessage(hwnd, idc, CB_SETCURSEL, seln, 0);
}

/**
 * Returns 1 if palette is one of the builtins, else 0
 */
int get_builtin_palette(void)
{
	int tmp, sz = NUM_ELEM(palette_strs);
	char str[256];
	man_calc_struct* m = &main_man_calc_struct;

	GetDlgItemText(hwnd_dialog, IDC_PALETTE, str, sizeof(str));

	if ((tmp = get_string_index(str, palette_strs, sz)) >= 0)
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

	GetDlgItemText(hwnd_dialog, IDC_PALETTE, pal_file, sizeof(pal_file));

	// See whether it's a BMP or text palette file; load using corresponding 
	// function. This is safe because files have to have 3-char extensions to 
	// make it into the dropdown list.
	bmp_flag = !_strnicmp(pal_file + strlen(pal_file) - 3, "bmp", 3);

	if ((fp = open_file(pal_file, "", bmp_flag)) != NULL)
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

int get_rendering(void)
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

int get_algorithm(void)
{
	char str[256];

	GetDlgItemText(hwnd_dialog, IDC_ALGORITHM, str, sizeof(str));
	return get_string_index(str, algorithm_strs, NUM_ELEM(algorithm_strs));
}

void get_num_threads(void)
{
	char str[256];
	man_calc_struct* m = &main_man_calc_struct;

	GetDlgItemText(hwnd_dialog, IDC_THREADS, str, sizeof(str));
	m->threads_i = get_string_index(str, threads_strs, NUM_ELEM(threads_strs));

	m->num_threads = 1 << m->threads_i;
	m->num_stripes = get_num_stripes(m->threads_i);
}

/**
 * Get all fields from the dialog.
 */
void get_dialog_fields(void)
{
	man_calc_struct* m = &main_man_calc_struct;

	m->max_iters = GetDlgItemInt(hwnd_dialog, IDC_ITERS, NULL, FALSE);

	update_iters(m, 0, 0); // clip
	SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);

	// Doing these only in the dialog box handler mean they don't take effect 
	// until user actually clicks them, whereas getting them here returns the 
	// realtime values.

	// It's useful to be able to set these on the fly
	m->rendering = get_rendering();
	m->precision = get_precision();
	m->algorithm = get_algorithm();

	// don't allow on-the-fly change to this (unimplemented)
	if (m->precision == PRECISION_EXTENDED)
		m->precision = PRECISION_DOUBLE;

	get_builtin_palette();

	// This not so much
	// get_num_threads();
}

/**
 * Update the slider position and get the new value. Clips any bad values.
 */
int set_slider_pos(int dlg_item, int pos)
{
	SendDlgItemMessage(hwnd_dialog, dlg_item, TBM_SETPOS, TRUE, (LONG)pos);
	return (int)SendDlgItemMessage(hwnd_dialog, dlg_item, TBM_GETPOS, 0, 0);
}

void setup_sliders(void)
{
	cfg.pan_rate.c_val  = set_slider_pos(IDC_PAN_RATE, cfg.pan_rate.c_val);
	cfg.zoom_rate.c_val = set_slider_pos(IDC_ZOOM_RATE, cfg.zoom_rate.c_val);
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
	if (!(m->status & STAT_DOING_SAVE))
	{
		sprintf_s(s, sizeof(s), "%s%s", calc ? "Calculating..." : "Ready ",
			calc ? "" : m->precision_loss ? "[Prec Loss]" : "");

		SetWindowText(hwnd_status_1, s);
	}

	sprintf_s(s, sizeof(s), "%d/%d  %c", log_pos + 1, log_count,
		m->precision == PRECISION_SINGLE
		? 'S' : m->precision == PRECISION_DOUBLE ? 'D' : 'E');

	SetWindowText(hwnd_status_2, s);
}

/**
 * If the palette has been modified in some way (inverted, locked, colors
 * changed with xor, etc) print a '*' before the Palette text.
 */
void print_palette_status()
{
	char s[32];
	man_calc_struct* m = &main_man_calc_struct;

	sprintf_s(s, sizeof(s), "%c Palette",
		(m->status & STAT_PALETTE_LOCKED) || m->pal_xor ? '*' : ' ');

	SetWindowText(GetDlgItem(hwnd_dialog, IDC_PAL_TEXT), s);
}

void not_implemented_yet(void)
{
	const char str[] = { "This feature is not implemented yet." };

	MessageBox(NULL, str, NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
}

void unsupported_alg_prec(void)
{
	const char str[] = {
		"Your (obsolete) CPU cannot run this algorithm/precision combination.\n"
		"Using C algorithm."
	};

	MessageBox(NULL, str, NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
}

int unrecommended_alg(void)
{
	const char str[] = {
		"Using the Fast algorithm with Normalized rendering may\n"
		"cause image artifacts. Switch to the Exact algorithm?"
	};

	return MessageBox(NULL, str, "Warning", 
		MB_YESNO | MB_ICONWARNING | MB_TASKMODAL);
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
	if (!(m->algorithm & ALG_EXACT) && (m->rendering == REN_NORMALIZED))
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

	m->rendering = get_rendering();
	m->precision = get_precision();
	m->algorithm = get_algorithm();

	if (m->precision == PRECISION_EXTENDED)
	{
		not_implemented_yet();
		SendDlgItemMessage(hwnd, IDC_PRECISION, CB_SETCURSEL,
			PRECISION_DOUBLE, 0);
	}
	else if (m->precision == PRECISION_DOUBLE)
	{
		// If CPU doesn't support SSE2, can only run C version
		if (m->sse_support < 2 && !(m->algorithm & ALG_C))
		{
			unsupported_alg_prec();
			SendDlgItemMessage(hwnd, IDC_ALGORITHM, CB_SETCURSEL,
				ALG_FAST_C, 0);
		}
	}
	else
	{
		// If CPU doesn't support SSE, can only run C version
		if (!m->sse_support && !(m->algorithm & ALG_C))
		{
			unsupported_alg_prec();
			SendDlgItemMessage(hwnd, IDC_ALGORITHM, CB_SETCURSEL,
				ALG_FAST_C, 0);
		}
	}

	// Give user a choice to switch to exact alg if using normalized rendering

	if (!(m->algorithm & ALG_EXACT) && (m->rendering == REN_NORMALIZED))
		if (unrecommended_alg() == IDYES)
		{
			SendDlgItemMessage
			(hwnd, IDC_ALGORITHM, CB_SETCURSEL, m->algorithm |= ALG_EXACT, 0);

			// need to recalc if switching to exact
			m->status |= STAT_RECALC_FOR_PALETTE;
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
		m->status |= STAT_NEED_RECALC;
	if (m->status & STAT_NEED_RECALC)	// if need recalculation,
		recalc_all = 1;	// force info update and wait cursor
	if (recalc_all)
	{
		// could be recalculating for palette- renormalize
		update_re_im(m, m->pan_xoffs, m->pan_yoffs);

		reset_quadrants(m);	// reset quadrants to UL only
		if (!do_rtzoom)	// If not realtime zooming
		{
			// print status line and reset fps timing values
			print_status_line(1);
			reset_fps_values(m);
			cursor = GetCursor();	// save current cursor
			SetCursor(wait_cursor);	// set wait cursor
		}
		// do this always, so we can change max_iters while zooming, etc.
		get_dialog_fields();

		// no longer need to recalc for palette
		m->status &= ~STAT_RECALC_FOR_PALETTE;
	}

	man_calculate_quadrants(m);

	// last iters actually calculated, for palette code
	m->max_iters_last = m->max_iters;

	InvalidateRect(hwnd_main, NULL, 0); // cause repaint with image data
	UpdateWindow(hwnd_main);

	// Don't update this stuff if realtime zooming. Will be done at intervals
	if (recalc_all && !do_rtzoom)
	{
		print_image_info(1);	// update image info
		print_status_line(0);

		SetCursor(cursor);	// restore old cursor
	}
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
			ypos += y_border;
		}
	}
	else // keep at current position
	{
		xpos = rc_dialog.left;
		ypos = rc_dialog.top;
	}

	SetWindowPos(hwnd_dialog, HWND_TOP, xpos, ypos, 0, 0, 
		SWP_NOSIZE | (hide ? SWP_HIDEWINDOW : SWP_SHOWWINDOW));
}

/**
 * Need to calculate all this stuff to accomodate strange window borders, large
 * fonts and such.
 */
void get_system_metrics(void)
{
	x_border = 2 * GetSystemMetrics(SM_CXSIZEFRAME);
	y_border = 2 * GetSystemMetrics(SM_CYSIZEFRAME) +
		GetSystemMetrics(SM_CYCAPTION);
}

/**
 * Enter and exit fullscreen mode.
 */
void toggle_fullscreen()
{
	man_calc_struct* m = &main_man_calc_struct;

	if ((m->status ^= STAT_FULLSCREEN) & STAT_FULLSCREEN)
	{
		// Hide dialog on entry to fullscreen mode, if so configured
		if (!(cfg.options.c_val & OPT_FULLSCREEN))
			m->status |= STAT_DIALOG_HIDDEN;

		SetWindowLongPtr(hwnd_main, GWL_STYLE, WS_POPUP | WS_VISIBLE);
		SetWindowPos(hwnd_main, NULL, 0, 0, 
			GetSystemMetrics(SM_CXSCREEN),
			GetSystemMetrics(SM_CYSCREEN), 
			SWP_DRAWFRAME | SWP_NOZORDER);
	}
	else // exit
	{
		// always show dialog when coming out of fullscreen
		m->status &= ~STAT_DIALOG_HIDDEN;

		SetWindowLongPtr(hwnd_main, GWL_STYLE,
			WS_OVERLAPPEDWINDOW | WS_VISIBLE);
		SetWindowPos(hwnd_main, NULL, main_rect.left, main_rect.top,
			main_rect.right - main_rect.left,
			main_rect.bottom - main_rect.top,
			SWP_NOZORDER);
	}
	UpdateWindow(hwnd_main);

	// recalc after resize, if enabled
	if (cfg.options.c_val & OPT_RECALC_ON_RESIZE)
		m->status |= STAT_RECALC_IMMEDIATELY;
}

/**
 * Resize window if necessary. Don't allow window size changes in fullscreen
 * mode, except by restore command (assume if user went to fullscreen, he wants
 * to stay there).
 *
 * v1.07
 */
void resize_window()
{
	man_calc_struct* m = &main_man_calc_struct;

	if (!(m->status & STAT_FULLSCREEN))
	{
		// xsize 0: fullscreen (other values < MIN_SIZE can be used later)
		if (cfg.xsize.c_val < MIN_SIZE)
			toggle_fullscreen();
		else if (cfg.ysize.c_val >= MIN_SIZE && // ignore any restoring ysize = 0
				(cfg.xsize.c_val != prev_xsize ||
				 cfg.ysize.c_val != prev_ysize))
		{
			// Without this, you can change the size of a maximized window, 
			// which puts it into a bad state (can't resize)
			ShowWindow(hwnd_main, SW_RESTORE);
			SetWindowPos(hwnd_main, HWND_TOP, 0, 0,
				cfg.xsize.c_val + x_border,
				cfg.ysize.c_val + y_border, 
				SWP_NOMOVE | SWP_NOCOPYBITS);

			// not sure why this is necessary. Sizes get out of sync?
			UpdateWindow(hwnd_main);

			prev_xsize = cfg.xsize.c_val;	// Save previous size
			prev_ysize = cfg.ysize.c_val;	// Can't do this in create_bitmap
		}
	}
	else if (cfg.ysize.c_val < MIN_SIZE)
		toggle_fullscreen(); // If in fullscreen mode: only allow restore
}

void png_error_msg(void)
{
	MessageBox(NULL,
		"The PNG library returned an error. Possible causes:\n\n"
		"1. The image you are trying to save may be too large.\n"
		"2. The filename for the image is invalid.",
		NULL, MB_OK | MB_ICONSTOP | MB_TASKMODAL);
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

	double start_time, t;
	man_calc_struct* m, * s;

	m = &main_man_calc_struct;
	s = &save_man_calc_struct;

	// Get sizes for image save, and clip
	if ((save_xsize = m->xsize) < MIN_SIZE)
		save_xsize = MIN_SIZE;
	if ((save_ysize = m->ysize) < MIN_SIZE)
		save_ysize = MIN_SIZE;

	// Get filename and add .png extension if not already present
	GetDlgItemText(hwnd_dialog, IDC_IMAGEFILE, img_file, sizeof(img_file));
	n = (int)strlen(img_file);
	if (n < 4 || _strnicmp(&img_file[n - 4], ".png", 4))
		strcat_s(img_file, sizeof(img_file), ".png");

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

	// Exact will be faster for 1-pixel high rows. 
	// Want for best quality anyway.
	s->algorithm = m->algorithm | ALG_EXACT;

	s->palette = m->palette;
	s->prev_pal = 0xFFFFFFFF;			// always recalc. pal lookup table 
										// before starting
	s->pal_xor = m->pal_xor;
	s->max_iters_color = m->max_iters_color;
	s->rendering = m->rendering;
	s->flags |= FLAG_CALC_RE_ARRAY;		// tell man_calculate to calculate the 
										// real array initially
	s->sse_support = m->sse_support;
	s->num_stripes = m->num_stripes;
	s->num_threads = m->num_threads;

	// Make sure all image data above is already captured before possibly 
	// popping a message box. Because this func is in a separate thread, the 
	// image can be modified (panned, etc) while the box is open.

	// If file already exists, ask user to confirm overwite
	if (!fopen_s(&fp, img_file, "rb"))
	{
		fclose(fp);
		sprintf_s(c, sizeof(c), "%s already exists. Overwrite?", img_file);
		if (MessageBox(NULL, c, "Warning",
			MB_YESNO | MB_ICONWARNING | MB_TASKMODAL) != IDYES)
		{
			m->status &= ~STAT_DOING_SAVE;
			return 0;
		}
	}

	if (!png_save_start(img_file, save_xsize, save_ysize))
	{
		m->status &= ~STAT_DOING_SAVE;
		png_error_msg();
		return 0;
	}

	free_man_mem(s); // free any existing arrays and alloc new arrays
	alloc_man_mem(s, save_xsize, 1);

	start_time = t = get_timer();

	s->pan_yoffs = -((save_ysize - 1) >> 1); // pan_yoffs of image top
	img_save = 1;

	for (i = 0; i < save_ysize; i++)
	{
		// iterate the row
		s->man_calculate(0, save_xsize - 1, 0, 0);

		// don't have to recalculate real array on subsequent rows
		s->flags &= ~FLAG_CALC_RE_ARRAY;

		// Palette-map the iteration counts to RGB data in png_buffer. 
		// Magnitudes are also available here.
		apply_palette((unsigned int*)s->png_buffer, s->iter_data,
			save_xsize, 1);

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
		{
			png_error_msg();
			break;
		}

		s->pan_yoffs++; // go to next row

		// print progress indicator every 0.5s
		if (get_seconds_elapsed(t) > 0.5)
		{
			sprintf_s(c, sizeof(c), "Saving... (%3.1f%%)",
				100.0 * (double)i / (double)save_ysize);
			SetWindowText(hwnd_status_1, c);
			t = get_timer();
		}
	}

	if (!png_save_end())
		png_error_msg();

	sprintf_s(c, sizeof(c), "Saved in %.1fs", get_seconds_elapsed(start_time));
	SetWindowText(hwnd_status_1, c);

	m->status &= ~STAT_DOING_SAVE;
	img_save = 0;
	return 1;
}

/**
 * The help message box.
 */
unsigned int CALLBACK show_help(LPVOID param)
{
	man_calc_struct* m = &main_man_calc_struct;

	const char str[] = {
		"For complete documentation, please go to the QuickMAN project webpage\n"
		"and click on the Documentation tab.\n\n"

		"https://sourceforge.net/projects/quickman/\n\n"

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

	// Using MB_TOPMOST gives an old window style... need to use SYSTEMMODAL, 
	// but then we get an annoying icon on the title bar.
	MessageBox(NULL, str, "QuickMAN Help", MB_OK | MB_SYSTEMMODAL);

	m->status &= ~STAT_HELP_SHOWING;
	return 0;
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

//#define STRIPE_THICKNESS	16 // this now comes from global config setting

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

	this_y = cfg.blit_stripe_thickness.c_val; // STRIPE_THICKNESS;
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

/**
 * Microsoft code for confining the mouse cursor to the main window.
 */
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

void fancy_intro()
{
#define MAG_STEP	1.07 // slow down a bit from previous versions

	man_calc_struct* m = &main_man_calc_struct;

	set_home_image(m);
	SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
	m->max_iters = 64;
	m->mag = MAG_START;
	do_rtzoom = 1; // prevent status line from being updated

	do
	{
		do_man_calculate(1);
		m->mag *= MAG_STEP;
	} while (m->mag <= 1.35);

	set_home_image(m);
	SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
	do_man_calculate(1);
	do_rtzoom = 0;
	
	// "resized" initially, but no need to calc again
	m->status &= ~STAT_RECALC_IMMEDIATELY;

	print_image_info(1);	// update image info
	print_status_line(0);	// update status bar

	m->calc_time = 0.0; // don't count intro time in file total time
}

/**
 * Initialize values that never change. Call once at the beginning of the
 * program.
 */
void man_init(void)
{
	int i, j;
	void* e;
	man_pointstruct* ps_ptr;
	man_calc_struct* m;

	// Initialize the calculation structure
	for (j = 0; j < 2; j++)
	{
		m = j ? &save_man_calc_struct : &main_man_calc_struct;

		m->flags = j ? FLAG_IS_SAVE | FLAG_CALC_RE_ARRAY : FLAG_CALC_RE_ARRAY;
		m->palette = DEFAULT_PAL;
		m->rendering = cfg.options.c_val & OPT_NORMALIZED
			? REN_NORMALIZED
			: REN_STANDARD;
		m->precision = PRECISION_AUTO;
		m->mag = HOME_MAG;
		m->max_iters = HOME_MAX_ITERS;

		// Initialize the thread state structures
		for (i = 0; i < MAX_THREADS; i++)
		{
			m->thread_states[i].thread_num = i;
			m->thread_states[i].calc_struct = m;

			// Create an auto-reset done event for each thread. 
			// The thread sets it when done with a calculation.
			e = CreateEvent(NULL, FALSE, FALSE, NULL);
			m->thread_states[i].done_event = e;
			m->thread_done_events[i] = e;

			// Init each thread's point structure
			m->thread_states[i].ps_ptr = ps_ptr = &m->pointstruct_array[i];

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

/**
 * Handler for the control dialog box. How do we make this stop blocking the
 * rest of the code when the user is moving it?
 *
 * sdiag
 */
INT_PTR CALLBACK
man_dialog_proc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	static int adj_iters_prev = 0;
	static int ignore_next_change = 0;
	static int new_file_entered = 0;
	static int new_file_selected = 0;
	int tab_spacing = 26; // dialog box units
	log_entry* e;
	man_calc_struct* m = &main_man_calc_struct;

	switch (uMsg)
	{
	case WM_INITDIALOG:
	{
		// Initialize all dialog box fields
		SetDlgItemInt(hwnd, IDC_ITERS, m->max_iters, FALSE);
		SendDlgItemMessage(hwnd,
			IDC_PAN_RATE,
			TBM_SETRANGE, TRUE,
			MAKELONG(0, MAX_PAN_RATE));

		SendDlgItemMessage(hwnd,
			IDC_ZOOM_RATE,
			TBM_SETRANGE, TRUE,
			MAKELONG(0, MAX_ZOOM_RATE));

		combobox_init(hwnd, IDC_PALETTE);
		combobox_init(hwnd, IDC_RENDERING);
		combobox_init(hwnd, IDC_PRECISION);
		combobox_init(hwnd, IDC_ALGORITHM);
		combobox_init(hwnd, IDC_THREADS);
		combobox_init(hwnd, IDC_LOGFILE);

		// default save filename- later have this scan for next available
		SetWindowText(GetDlgItem(hwnd, IDC_IMAGEFILE), img_file);

		// Set tab stops for the INFO window
		SendDlgItemMessage(hwnd, IDC_INFO, EM_SETTABSTOPS, 1,
			(LPARAM)&tab_spacing);

		// If any of these are null, the system has major problems
		hwnd_iters		= GetDlgItem(hwnd, IDC_ITERS);
		hwnd_info		= GetDlgItem(hwnd, IDC_INFO);
		hwnd_status_1	= GetDlgItem(hwnd, IDC_STATUS_1);
		hwnd_status_2	= GetDlgItem(hwnd, IDC_STATUS_2);

		return FALSE;
	}
	case WM_VSCROLL: // Update iterations from spin control
	{
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
				{
					update_iters(m, 0, 1);
					SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
				}
				else
				{
					update_iters(m, 1, 0);
					SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
				}
				adj_iters_prev = HIWORD(wParam);
				break;
			}
		}
		return TRUE;
	}
	case WM_HSCROLL:
	{
		// Slider values are being updated. Moving these around during zooming 
		// annoyingly creates about a 1% slowdown, even with no processing.

		if ((HWND)lParam == GetDlgItem(hwnd, IDC_PAN_RATE))
			cfg.pan_rate.c_val = (int)SendDlgItemMessage
			(hwnd, IDC_PAN_RATE, TBM_GETPOS, 0, 0);

		if ((HWND)lParam == GetDlgItem(hwnd, IDC_ZOOM_RATE))
			cfg.zoom_rate.c_val = (int)SendDlgItemMessage
			(hwnd, IDC_ZOOM_RATE, TBM_GETPOS, 0, 0);

		// reset frames/sec timing values when pan or zoom rate changes
		reset_fps_values(m);

		return TRUE;
	}
	case WM_COMMAND:
	{
		switch (LOWORD(wParam))
		{
		case IDC_LOGFILE:
		{
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
		}
		case IDC_PALETTE:
		case IDC_RENDERING:
		{
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
				if ((m->status & STAT_RECALC_FOR_PALETTE) ||
					(m->max_iters != m->max_iters_last))
				{
					// update re/im from any pan offsets and reset offsets
					update_re_im(m, m->pan_xoffs, m->pan_yoffs);
					do_man_calculate(1);
				}

				// Apply palette to the whole image (in UL quadrant here)
				apply_palette(m->quad[UL].bitmap_data,
					m->iter_data,
					m->xsize,
					m->ysize);

				InvalidateRect(hwnd_main, NULL, 0); // cause repaint
				UpdateWindow(hwnd_main);
			}
			return TRUE;
		}
		case IDC_ALGORITHM:
		case IDC_PRECISION:
		{
			if (HIWORD(wParam) == CBN_SELCHANGE)
				// make sure CPU supports alg/precision combination
				check_alg(hwnd);
			return TRUE;
		}
		case IDC_THREADS:
		{
			if (HIWORD(wParam) == CBN_SELCHANGE)
				get_num_threads();
			return TRUE;
		}
		case ID_HOME:
		{
			// Reset to base image coordinates; deliberate fallthru
			set_home_image(m);
			SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);

			// only do this on home or log next/prev, not recalculation
			autoreset_settings(&cfg);

			resize_window();
		}
		case ID_CALCULATE:
		{
			// stop any rt zoom in progress
			do_rtzoom = prev_do_rtzoom = 0;

			// update re/im from any pan offsets and reset offsets
			update_re_im(m, m->pan_xoffs, m->pan_yoffs);
			reset_pan_state(m);		// reset pan filters and movement state
			get_pan_steps(NULL, NULL, 0);	// reset any pan lock
			print_palette_status();
			do_man_calculate(1);	// calculate all pixels, update image info
			SetFocus(hwnd_main);	// allow arrow keys to work immediately for 
									// panning
			// set focus AFTER calculating, or you get an annoying blink
			return TRUE;
		}
		case ID_LOG_IMAGE:
		case ID_LOG_PREV:
		case ID_LOG_NEXT:
		{
			// Get the current filename
			GetDlgItemText(hwnd, IDC_LOGFILE, log_file, sizeof(log_file));

			// If it's new, read it. If it was a newly entered filename, add it 
			// to the list.
			if (new_file_entered || new_file_selected)
			{
				if (LOWORD(wParam) != ID_LOG_IMAGE)
				{
					log_read(log_file, "", 1);	// read logfile into array

					// for testing load-balancing alg
					reset_thread_load_counters(m);
				}
				if (new_file_entered)
					SendDlgItemMessage(hwnd, IDC_LOGFILE, CB_ADDSTRING, 0,
						(LPARAM)log_file);
			}
			if (LOWORD(wParam) == ID_LOG_IMAGE)
			{
				// Need to update re/im with panning offsets, or logged 
				// coordinates will be wrong!
				// Added this line for v1.03 bug fix
				update_re_im(m, m->pan_xoffs, m->pan_yoffs);

				// reset pos if new file
				log_update(log_file, new_file_entered | new_file_selected);

				// update current/total number of log images
				print_status_line(0);

				// Ok if this stays till next calculation
				SetWindowText(hwnd_status_1, "Logged");
			}
			new_file_entered = 0;
			new_file_selected = 0;

			// If user wants a new log image...
			if ((LOWORD(wParam) == ID_LOG_NEXT || LOWORD(wParam == ID_LOG_PREV)) 
				&& log_count)
			{
				// Autoreset any previous settings that need it
				autoreset_settings(&cfg);

				// Clear any pan offsets (no need to update re/im as they are 
				// reset from logfile)
				m->pan_xoffs = 0;
				m->pan_yoffs = 0;

				// get next/prev image from logfile array
				if ((e = log_get(LOWORD(wParam) == ID_LOG_NEXT)) == NULL)
					return TRUE;

				// Update any new settings: 0 = don't copy to default_val
				copy_changed_settings(&cfg, &e->log_settings, 0);

				// Update sliders, info box, iters, and palette
				setup_sliders();
				update_iters(m, 0, 0);
				SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
				UpdateWindow(hwnd_iters);

				// stop any rt zoom in progress
				do_rtzoom = prev_do_rtzoom = 0;

				reset_pan_state(m);	// reset pan filters and movement state

				// set pan lock (or turn it off if pan_key is 0)
				get_pan_steps(NULL, NULL, cfg.pan_key.c_val);

				// Change palette, max iters color, and inversion status if not 
				// locked
				if (!(m->status & STAT_PALETTE_LOCKED))
				{
					m->pal_xor = cfg.pal_xor.c_val;
					m->max_iters_color = cfg.max_iters_color.c_val;

					SendDlgItemMessage(hwnd, IDC_PALETTE, CB_SETCURSEL,
						m->palette, 0);

					// if user palette, read it in
					if (m->palette >= num_palettes)
						get_user_palette();
				}

				print_palette_status();
				resize_window();

				// If not zooming, just calculate
				if (!cfg.zoom_in_out.c_val)
				{
					// calculate all pixels, update image info
					do_man_calculate(1);

					// already calculated: don't need to do again
					m->status &= ~STAT_RECALC_IMMEDIATELY;

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
		}
		case ID_ZOOM:	// Do a realtime zoom in to the current image
		{
			// update re/im from any pan offsets and reset offsets
			update_re_im(m, m->pan_xoffs, m->pan_yoffs);

			reset_fps_values(m);
			reset_thread_load_counters(m);
			zoom_start_time = get_timer();
			zoom_start_mag = m->mag;
			m->mag = MAG_MIN;
			do_rtzoom = RTZOOM_IN | RTZOOM_BTN;
			return TRUE;
		}
		case ID_FULLSCREEN:
		{
			toggle_fullscreen();
			return TRUE;

			// case ID_SLIDESHOW:
			// case ID_OPTIONS:
			   // not_implemented_yet();
			   // return TRUE;
		}
		case ID_SAVE_IMAGE:
		{
			if (!(m->status & STAT_DOING_SAVE))
			{
				// prevent re-entering function when already saving
				m->status |= STAT_DOING_SAVE;

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
		}
		case ID_HELP_BUTTON: // MS won't generate ID_HELP
		{
			// Have to do this in a separate thread just to keep the messagebox 
			// from blocking my main window...

			// don't show help box if it's already showing
			if (!(m->status & STAT_HELP_SHOWING))
			{
				m->status |= STAT_HELP_SHOWING;
				QueueUserWorkItem(show_help, NULL,
					WT_EXECUTELONGFUNCTION | (MAX_QUEUE_THREADS << 16));
			}
			return TRUE;
		}
		default:
			return FALSE;
		}
		break;
	}
	case WM_CLOSE:	// Never close or destroy this dialog
		break;
	case WM_DESTROY:
		break;
	}
	return FALSE;	// return FALSE if we didn't process the message 
					// (exceptions for INITDIALOG) otherwise TRUE
}

/**
 * The window function for the main window.
 */
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
	static int prev_nav_mode = NAV_RTZOOM;
	static int prev_sizing = 0;
	static int prev_max_restore = 0;	// used to detect when window 
										// transitioned from maximized to 
										// restored

	man_calc_struct* m = &main_man_calc_struct;

	switch (nMsg)
	{
	case WM_CREATE: // The window is being created
	{
		hscreen_dc = CreateCompatibleDC(NULL); // screen device context
		hwnd_dialog = CreateDialog
		(hinstance, MAKEINTRESOURCE(IDD_MAN_DIALOG), hwnd, man_dialog_proc);

		// pen for the zoom rectangle; easier to see than PS_DOT
		hpen = CreatePen(PS_SOLID, 2, RGB(0, 0, 0));

		setup_sliders();
		set_alg_warning(); // warn if normalized rendering and fast alg

		return FALSE;
	}
	case WM_PAINT:	// The window needs to be painted (redrawn).
	{
		hdc = BeginPaint(hwnd, &ps);

		// Blit the mandelbrot bitmap. Could take rectangular regions from 1 to 
		// 4 quadrants.
		// Blits thin horizontal stripes. See comments at striped_blit().
		// Do upper quadrants (UL and UR), then lower (LL and LR).

		// Also optimize for the case where we've recalculated the whole image 
		// (only UL is valid). In this case do a normal blit.
		q = &m->quad[UL];
		if (q->blit_xsize == m->xsize && q->blit_ysize == m->ysize)
		{
			SelectObject(hscreen_dc, q->handle);
			BitBlt(hdc, 0, 0, 
				q->blit_xsize, 
				q->blit_ysize, 
				hscreen_dc, 0, 0, SRCCOPY);
		}
		else
		{
          	striped_blit(&m->quad[UL], &m->quad[UR], hdc, hscreen_dc);
			striped_blit(&m->quad[LL], &m->quad[LR], hdc, hscreen_dc);
		}
		EndPaint(hwnd, &ps);
		return FALSE;
	}
	case WM_SIZING:
	{
		m->status |= STAT_RECALC_IMMEDIATELY;
		return TRUE;
	}
	case WM_LBUTTONDOWN:
	{
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

		if (nav_mode == NAV_PAN)
			SetCursor(hclosed_cursor);	// set closed hand
		if (nav_mode == NAV_RTZOOM)	// if in realtime zoom mode,
			do_rtzoom = RTZOOM_IN;	// set global flag for do_zooming()
		else
			allow_mode_change = 0;		// don't allow mode change while button 
										// is down if not in rtzoom mode

		return FALSE;
	}
	case WM_LBUTTONUP: // Zoom in
	{
		mouse_x[1] = LOWORD(lParam);
		mouse_y[1] = HIWORD(lParam);
		allow_mode_change = 1; // Allow mode change again after button released

		if (dragging)
		{
			dragging = 0;
			if (nav_mode == NAV_ZOOM)
			{
				// Update mandelbrot parms from rectangle and recalculate
				update_re_im_mag(m, have_box, 1,
					mouse_x[0],
					mouse_y[0],
					mouse_x[1],
					mouse_y[1]);

				do_man_calculate(1);
			}
			else
			{
				// Update image info (excluding iters/sec) after pan or 
				// realtime zoom
				print_image_info(0);
			}
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
		if (nav_mode != NAV_PAN)	// this implements zoom lock
			do_rtzoom = prev_do_rtzoom = 0;	// clear realtime zoom 
											// flag for do_zooming()

		ClipCursor(NULL); // release mouse cursor and capture
		ReleaseCapture();

		return FALSE;
	}
	case WM_MOUSEMOVE:
	{
		mouse_x[1] = LOWORD(lParam);
		mouse_y[1] = HIWORD(lParam);

		// Get re/im coords at the mouse position, for realtime zoom
		get_mouse_re_im(mouse_x[1], mouse_y[1]);

		// If user is dragging, draw the zoom rectangle (zoom mode) or pan the 
		// image (pan mode)
		if (nav_mode == NAV_PAN) // panning mode- move the image
		{
			// allow panning using right button drag also
			if (wParam & (MK_LBUTTON | MK_RBUTTON))
			{
				int offs_x, offs_y;

				// Get difference from previous mouse location; use as pan 
				// offset
				offs_x = mouse_x[1] - mouse_x[0];
				offs_y = mouse_y[1] - mouse_y[0];

				// update previous mouse location
				mouse_x[0] = mouse_x[1];
				mouse_y[0] = mouse_y[1];

				pan_image(m, offs_x, offs_y);  // do the pan
				do_man_calculate(0);
			}
		}
		else if (nav_mode == NAV_ZOOM)
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
						prev_mouse_x, 
						prev_mouse_y);

					have_box = 1;
				}

				prev_mouse_x = mouse_x[1];
				prev_mouse_y = mouse_y[1];

				// draw new rectangle
				Rectangle(hdc, mouse_x[0], mouse_y[0], mouse_x[1], mouse_y[1]);
				ReleaseDC(hwnd, hdc);
			}
		}
		return FALSE;
	}
	case WM_RBUTTONDOWN: // Zoom out
	{
		mouse_x[0] = LOWORD(lParam);
		mouse_y[0] = HIWORD(lParam);
		mouse_x[1] = LOWORD(lParam); // Init this for switch from zoom to pan also
		mouse_y[1] = HIWORD(lParam);

		// update re/im from any pan offsets and reset offsets
		update_re_im(m, m->pan_xoffs, m->pan_yoffs);

		// get re/im coords at the mouse position, for realtime zoom
		get_mouse_re_im(mouse_x[0], mouse_y[0]);

		if (nav_mode == NAV_RTZOOM)	// if in realtime zoom mode,
			do_rtzoom = RTZOOM_OUT;	// set global flag for do_zooming()
		else
			allow_mode_change = 0;		// else don't allow mode change while 
										// button is down

		if (nav_mode == NAV_PAN)
			SetCursor(hclosed_cursor);	// set closed hand

		confine_mouse_cursor();	// also need to confine here, for realtime zoom
		return FALSE;
	}
	case WM_RBUTTONUP: // Zoom out
	{
		mouse_x[1] = LOWORD(lParam);
		mouse_y[1] = HIWORD(lParam);
		allow_mode_change = 1; // allow mode change again after button released

		// Zoom out from current point, and recenter
		if (nav_mode == NAV_ZOOM)
		{
			update_re_im_mag(m, 0, 0,
				mouse_x[0],
				mouse_y[0],
				mouse_x[1],
				mouse_y[1]);

			do_man_calculate(1);
		}

		// this implements zoom lock
		if (nav_mode != NAV_PAN)
			do_rtzoom = prev_do_rtzoom = 0;	// clear realtime zoom 
											// flag for do_zooming()

		ClipCursor(NULL);	// release mouse cursor and capture
		ReleaseCapture();
		return FALSE;
	}
	case WM_MOUSEWHEEL:	// Use mousewheel to adjust iterations. (Maybe palette 
						// too, if button down?)
	{
		if (GET_WHEEL_DELTA_WPARAM(wParam) > 0)
		{
			update_iters(m, 1, 0);
			SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
		}
		else
		{
			update_iters(m, 0, 1);
			SetDlgItemInt(hwnd_dialog, IDC_ITERS, m->max_iters, FALSE);
		}
		return FALSE;
	}
	case WM_KEYDOWN:	// Go to pan mode while space is held down (if allowed)
	{
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
					? NAV_ZOOM : NAV_RTZOOM;

				// Change cursor only if released in client area (kluge)
				if (GetCursor() != arrow_cursor)
					SetCursor(mag_zoom_cursor);
			}

			nav_mode = (wParam == VK_SPACE) ? NAV_PAN : prev_nav_mode;

			if (nav_mode == NAV_PAN)
			{
				// Reset mouse position for pan
				mouse_x[0] = mouse_x[1];
				mouse_y[0] = mouse_y[1];

				// save realtime zoom state
				prev_do_rtzoom = do_rtzoom;
				do_rtzoom = 0;	// stop any realtime zoom

				// Fix toggling between hand and arrow in non-client area 
				// (kluge)
				if (GetCursor() == mag_zoom_cursor)
					SetCursor(prev_do_rtzoom ?
						hclosed_cursor : hopen_cursor);
			}
		}

		// Handle the various hotkeys
		switch (wParam)
		{
		case 'C':	// 'C' toggles the control dialog on and off. Allow to work 
					// in non-fullscreen mode too.

			update_dialog((m->status ^= STAT_DIALOG_HIDDEN) &
				STAT_DIALOG_HIDDEN, 0); // 0 = don't move
			break;
		case VK_ESCAPE:	// ESC exits out of fullscreen mode but does not enter 
						// it.
			// deliberate fallthrough if fullscreen
			if (!(m->status & STAT_FULLSCREEN))
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

			m->status ^= STAT_PALETTE_LOCKED;
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
	}
	case WM_HELP:	// F1 shows help
	{
		SendMessage(hwnd_dialog, WM_COMMAND, ID_HELP_BUTTON, 0);
		return TRUE;
	}
	case WM_KEYUP: // Go back to zoom mode if space released, if allowed.
	{
		if (nav_mode == NAV_PAN)
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
	}
	case WM_SETCURSOR:	// We get this message whenever the cursor moves in 
						// this window.
	{
		// Set cursor (hand or zoom) based on nav mode. Also set keyboard focus 
		// to this window for key detection. Eliminates need to click window 
		// first to set focus
		SetFocus(hwnd);

		// only set zoom/hand cursor in client area
		if (LOWORD(lParam) == HTCLIENT)
		{
			SetCursor(nav_mode == NAV_PAN ? hopen_cursor : mag_zoom_cursor);
			return TRUE;
		}
		break; // let system set cursor outside client area (resize, arrow, etc)

		// This message comes after the user finishes a drag or movement, but, 
		// annoyingly, not after a maximize or restore. There's some extra code 
		// to handle those in WM_WINDOWPOSCHANGED
	}
	case WM_EXITSIZEMOVE:
	{
		// if previous resize/move was a resize, 
		if (prev_sizing)
			// recalculate image if enabled
			if (cfg.options.c_val & OPT_RECALC_ON_RESIZE)
			{
				// If we recalculate here, a partial bad-state window appears 
				// while the calculation is going on. Recalculate in main loop 
				// instead.
				m->status |= STAT_RECALC_IMMEDIATELY;
			}
		return FALSE;

		// Called on window sizing or changing position. Pretty wasteful to 
		// call create_bitmap (which frees and reallocates all arrays with 
		// every pixel of movement) here. Memory fragmentation? Changing this 
		// currently causes some issues. Maybe fix later.

		// Seems like we always get this at startup before the paint message.
	}
	case WM_WINDOWPOSCHANGED:
	{
		GetWindowPlacement(hwnd, &wp);

		// only do this if not minimizing window
		if (wp.showCmd != SW_SHOWMINIMIZED)
		{
			// Move dialog box along with main window
			update_dialog(m->status & STAT_DIALOG_HIDDEN, 1); // 1 = move

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
				print_image_info(0);	// update size info if resized
			}

			// Send a WM_EXITSIZEMOVE message if window was maximized/restored
			if (wp.showCmd != prev_max_restore)
			{
				SendMessage(hwnd, WM_EXITSIZEMOVE, 0, 0);
				prev_max_restore = wp.showCmd;
			}

			// Save main window rect for fullscreen mode if not maximized and 
			// not already in fullscreen
			if (wp.showCmd != SW_SHOWMAXIMIZED && !(m->status & STAT_FULLSCREEN))
				GetWindowRect(hwnd, &main_rect);

			return FALSE;
		}
		return TRUE;
	}
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

int CALLBACK
WinMain(HINSTANCE hInst, HINSTANCE hPrev, LPSTR lpCmd, int nShow)
{
	int i;
	MSG msg;
	WNDCLASSEX wndclass;
	static char* classname = "ManWin";
	man_calc_struct* m = &main_man_calc_struct;
	man_calc_struct* s = &save_man_calc_struct;

	m->man_calculate = &man_calculate;
	s->man_calculate = &man_calculate;

	m->apply_palette = &apply_palette;
	s->apply_palette = &apply_palette;

	hinstance = hInst;

	// Read any default settings from quickman.cfg. Do this early, because it 
	// can change almost anything used below.
	read_cfg_file();
	get_cpu_info();
	get_system_metrics();
	man_init();

	for (i = 1; i < MAX_THREADS; i++) // 0 is the master thread
	{
		m->pal_events[i] = CreateEvent(NULL, FALSE, FALSE, NULL);
		s->pal_events[i] = CreateEvent(NULL, FALSE, FALSE, NULL);

		if (m->pal_events[i] == NULL || s->pal_events[i] == NULL)
			return 0;
	}

	if (!(num_palettes = init_palettes(m, DIVERGED_THRESH)))
		return 0;

	// create a window class for our main window
	memset(&wndclass, 0, sizeof(WNDCLASSEX));
	wndclass.lpszClassName	= classname;
	wndclass.cbSize			= sizeof(WNDCLASSEX);
	wndclass.style			= CS_HREDRAW | CS_VREDRAW;
	wndclass.lpfnWndProc	= MainWndProc;
	wndclass.hInstance		= hInst;
	
	// use mini mandelbrot icon
	wndclass.hIcon			= LoadIcon(hInst, MAKEINTRESOURCE(IDI_MAN));

	// derive small icon from normal one we set other cursors based on nav mode
	wndclass.hIconSm		= NULL;

	wndclass.hCursor		= LoadCursor(NULL, IDC_ARROW);
	arrow_cursor			= wndclass.hCursor;
	wndclass.hbrBackground	= NULL;	// don't need this
	RegisterClassEx(&wndclass);

	// Load cursors - wait, zoom, hand open, and hand closed
	wait_cursor		= LoadCursor(NULL, IDC_WAIT);
	mag_cursor		= LoadCursor(hInst, MAKEINTRESOURCE(IDC_MAG));
	rtzoom_cursor	= LoadCursor(hInst, MAKEINTRESOURCE(IDC_RTZOOM));
	mag_zoom_cursor = rtzoom_cursor;
	hopen_cursor	= LoadCursor(hInst, MAKEINTRESOURCE(IDC_HAND_OPEN));
	hclosed_cursor	= LoadCursor(hInst, MAKEINTRESOURCE(IDC_HAND_CLOSED));

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

	UpdateWindow(hwnd_main);
	update_dialog(1, 1);          // 1 = hide, 1 = move
	ShowWindow(hwnd_main, nShow); // causes dialog to become visible
	UpdateWindow(hwnd_dialog);

	// Add user palettes and logfiles to their dropdown menus
	add_user_palettes_and_logfiles();

	// Read logfile - add_user_palettes() must be called before this
	log_read(LOG_FILE,
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
			{
				// don't use 100% of CPU when idle. Also see do_panning()
				Sleep(2);
			}
		}
	}

	free_man_mem(&main_man_calc_struct);
	free_man_mem(&save_man_calc_struct);

	return (int)msg.wParam;
}