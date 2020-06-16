/******************************************************************************
 * quickman.h
 * -> Header file for the QuickMAN SSE/SSE2-based Mandelbrot Set calculator
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
 * 11/12/08 PG: Initial for v1.10. Split off most definitions from the C files
 * into this header file.
 */

#ifndef _QUICKMAN_H_
#define _QUICKMAN_H_

#define CFG_FILE "quickman.cfg"	// configuration file containing default settings
#define LOG_FILE "quickman.log"	// default log file
#define IMG_FILE "image_1"		// default image file name
#define CFG_SIZE			13	// total number settings from configuration

/**
 * Precision used in calculation
 */
#define PRECISION_AUTO		0	// Automatically determined based on 
								// magnification
#define PRECISION_SINGLE	1	// 32-bit float
#define PRECISION_DOUBLE	2	// 64-bit double
#define PRECISION_EXTENDED	3	// 80-bit (x87) double or larger

/**
 * Available algorithms
 */
#define ALG_FAST_ASM_AMD	0	// Use the "wave" algorithm to guess pixels
#define ALG_EXACT_ASM_AMD	1	// Calculate every pixel (no interpolation or 
								// guessing)
#define ALG_FAST_ASM_INTEL	2	// Intel versions
#define ALG_EXACT_ASM_INTEL	3	//
#define ALG_FAST_C			4	// Unoptimized C versions
#define ALG_EXACT_C			5	//
//#define ALG_ERROR           6 // Show error image: exact ^ fast

#define ALG_EXACT			1	// using Exact alg if this bit set (change with 
								// above)
#define ALG_INTEL			2	// using Intel alg if this bit set
#define ALG_C				4	// using C alg if this bit set

/**
 * Rendering algorithms
 */
#define REN_STANDARD		0	// rendering algorithm standard (keep this 0)
#define REN_NORMALIZED		1	// rendering algorithm normalized

#define NUM_ELEM(a)			(sizeof(a) / sizeof(a[0]))

#define MAX_THREADS_IND		5	// Set this to set maximum number of threads 
								// (== 2^this)
#define MAX_THREADS			(1 << MAX_THREADS_IND)

/**
 * Max threads that can be running at once, with save going on in the 
 * background.
 */
#define MAX_QUEUE_THREADS  (MAX_THREADS * 2 + 3)

//#define USE_PERFORMANCE_COUNTER   // See get_timer()

#ifdef _WIN32
#include <windows.h>
#endif

#include <stdio.h>

#ifdef USE_PERFORMANCE_COUNTER
typedef LARGE_INTEGER TIME_UNIT;
#else
typedef DWORD TIME_UNIT;
#endif

#define DEFAULT_PAL		2		// new loud palette

/**
 * Home image parameters
 */
#define MAG_START		0.3
#define HOME_RE			-0.7
#define HOME_IM			0.001	// offset y a bit so x axis isn't black on the 
								// left side (looks a little better)
#define HOME_MAG		1.35
#define HOME_MAX_ITERS	256

/**
 * Navigation modes
 */
#define NAV_ZOOM	1	// zoom in/out 2x or zoom on rectangle, using magnifier
#define NAV_RTZOOM	2	// realtime zoom in/out
#define NAV_PAN		3	// pan around the image

/**
 * Bits for do_rtzoom
 */
#define RTZOOM_IN	1	// zoom in
#define RTZOOM_OUT	2	// zoom out
#define RTZOOM_BTN	4	// 1 if zooming with the zoom button, 
						// 0 if with the mouse

/**
 * Multiply or divide the magnification by this factor when the user
 * zooms using a single click.
 */
#define MAG_ZOOM_FACTOR		2.0
#define MAG_MIN				0.02	// minimum magnification

/**
 * Magnitude squared must exceed this for a point to be considered diverged.
 * Values smaller than 16 usually give less aesthetically pleasing results. 
 * Maybe make this user-settable.
 */

#define DIVERGED_THRESH		16.0	// 4.0
#define DIVERGED_THRESH_SQ	256		// DIVERGED_THRESH squared (integer)

/**
 * Integer portion of high word in double, corresponding to threshold above
 */
#define DIV_EXP				0x40300000	// for doubles with 11 bit exponent
#define DIV_EXP_FLOAT		0x41800000	// for floats with 8 bit exponent
//#define DIV_EXP			0x40100000	// for 4.0

#define MIN_ITERS	2			// allow to go down to min possible, for 
								// overhead testing
#define MAX_ITERS	0x08000000	// keep upper 4 bits free in iter array, in 
								// case we need them for something

#define MIN_SIZE	4			// min image size dimension. Code should work 
								// down to 1 x 1

//#define WM_DUMMY	(WM_APP + 100)	// use this form if messages needed

/**
 * Status bits
 */
#define STAT_NEED_RECALC		1   // 1 if image must be recalculated next 
									// time, for whatever reason
#define STAT_RECALC_FOR_PALETTE	2	// 1 if image must be recalculated before a 
									// new palette can be applied (true after 
									// panning)
#define STAT_FULLSCREEN			4	// 1 if in fullscreen mode
#define STAT_RECALC_IMMEDIATELY	8	// 1 if the image should be recalculated 
									// immediately (e.g., after window was 
									// resized)
#define STAT_DIALOG_HIDDEN		16	// 1 if the control dialog is currently 
									// hidden
#define STAT_PALETTE_LOCKED		32	// 1 if the palette is currently locked 
									// (ignore logfile palettes)
#define STAT_HELP_SHOWING		64	// 1 if the help window is showing
#define STAT_DOING_SAVE			128	// 1 if doing a save

/* ------------------- Quadrant-based panning structures ------------------- */

/**
 * Structure for rectangles used in panning
 */
typedef struct {
	int x[2];	// x coordinates; x[0] <= x[1]
	int y[2];	// y coordinates; y[0] <= y[1]
	int valid;	// nonzero if rect is valid, 0 if invalid (should be ignored)
} rectangle;

/**
 * Structure for quadrants used in panning
 */
typedef struct {
	int status;				// see below
	void* handle;			// handle to the bitmap

	rectangle quad_rect;	// rectangle coordinates of this quadrant

	// Raw bitmap data. Each 32-bit value is an RGB triplet (in bits 23-0) and 
	// one 0 byte (in bits 31-24). Faster to access than a 24-bit bitmap.
	unsigned int* bitmap_data;

	// Blitting parameters. All offsets are quadrant-relative (i.e., range from 
	// 0 to xsize - 1 and 0 to ysize - 1 inclusive).

	// Source offsets: blitting starts from this point in the bitmap
	int src_xoffs;
	int src_yoffs;

	// Dest offsets: blitting goes to this point on the screen
	int dest_xoffs;
	int dest_yoffs;
	
	// Size of area to blit
	int blit_xsize;
	int blit_ysize;
} quadrant;

/**
 * Values for the quadrant status word
 * If set, blit this quadrant's data to the screen
 */
#define QSTAT_DO_BLIT	1

#define UL 0	// upper left
#define UR 1	// upper right
#define LL 2	// lower left
#define LR 3	// lower right

/* -------------- Structures and variables used in iteration --------------- */

/**
 * Structure to hold the state of 4 iterating points (or 8 for SSE). For SSE 
 * (single precision), 32-bit floats are packed into the fields, otherwise 
 * 64-bit doubles. This needs to be 16-byte aligned. Project/compiler options 
 * can't guarantee this alignment: must use syntax below.
 *
 * Only the first 4/8 values of each array are used (unless otherwise noted). 
 * But the additional values are still necessary to force each array to occupy 
 * its own 64-byte cache line (i.e, no sharing). With line sharing there can be 
 * conflicts that cost cycles.
 *
 * May even want to give each 128 bits (xmm reg) its own cache line- change x 
 * to x01, x23, etc. Initialization and divergence detection would be nastier.
 *
 * Tried expanding x, y, and yy to 16 doubles so ..23 regs could have own cache 
 * line: no effect.
 */
typedef struct {
	double x[8];		// 0   x, y, yy = the state of the iterating points
	double y[8];		// 64
	double yy[8];		// 128
	double a[8];		// 192 Real coordinate of point
	double b[8];		// 256 Imag coordinate of point
	double mag[8];		// 320 Magnitudes from current iteration
	double magprev[8];	// 384 Magnitudes from previous iteration
	double two_d[8];	// 448 Only 1st 2 used; must be set to 2.0. Used in 
						// SSE2 routine.
	float two_f[16];	// 512 Only 1st 4 used; must be set to 2.0. Used in SSE 
						// routine.
	double rad_d[8];	// 576 Radius^2 for divergence detection; only 1st 2 
						// used; only used in Intel version
	float rad_f[16];	// 640 only 1st 4 used

	// Even though the following fields aren't used in the inner loop, there's 
	// a slight decrease in performance if they aren't aligned to a 64-byte 
	// cache line.

	// 704 Current iteration counts; only 1st 4/8 used
	unsigned int iters[16];
	
	// 768 Pointer into iteration count array; only 1st 4/8 used
	unsigned int* iters_ptr[16];
	
	// 832 Pointer into iteration count array; only 1st 4/8 used
	float* mag_ptr[16];
	
	// 896 Iterations counter, for benchmarking. M$ 64-bit aligns this, adding 
	// (crash-causing) extra padding if not already on a 64-bit boundary...
	unsigned long long iterctr;   
	
	// 904 loop sets ab_in to the point to iterate on 
	// (ab_in[0] = re, ab_in[1] = im). Others unused. MS also 64-bit aligns 
	// this
	double ab_in[2];
	
	unsigned int cur_max_iters;	// 920 Max iters to do this loop
	
	// Status of pointstruct queue (free/full slots)
	unsigned int queue_status;
	
	// Pad to make size a multiple of 64. Necessary, otherwise code will crash 
	// with 2 or more threads due to array misalignment.
	unsigned int pad[8];
} man_pointstruct;

/**
 * Pointers to structure members for asm functions.
 * There should be some function to calculate these automatically, but use 
 * constants for now. Changing structure order can slow things down anyway so 
 * it shouldn't be done without good reason.
 *
 * Aliases for 4 double-precision points. 
 * EBX points to the beginning of the structure.
 */
#define PS4_X01            [ebx + 0]
#define PS4_X23            [ebx + 0 + 16]
#define PS4_Y01            [ebx + 64]
#define PS4_Y23            [ebx + 64 + 16]
#define PS4_YY01           [ebx + 128]
#define PS4_YY23           [ebx + 128 + 16]
#define PS4_A01            [ebx + 192]
#define PS4_A23            [ebx + 192 + 16]
#define PS4_B01            [ebx + 256]
#define PS4_B23            [ebx + 256 + 16]
#define PS4_MAG01          [ebx + 320]		// Magnitudes of points 0 and 1

// Locations of exponent bits in magnitudes
#define PS4_MEXP0          [ebx + 320 + 4]

#define PS4_MEXP1          [ebx + 320 + 12]
#define PS4_MAG23          [ebx + 320 + 16]	// Magnitudes of points 2 and 3
#define PS4_MEXP2          [ebx + 320 + 20]
#define PS4_MEXP3          [ebx + 320 + 28]

// Magnitudes of points 0 and 1 after the previous iteration
#define PS4_MAGPREV01      [ebx + 384]

#define PS4_MAGPREV23      [ebx + 384 + 16]	// ditto for points 2 and 3
#define PS4_TWO            [ebx + 448]
#define PS4_RAD            [ebx + 576]
#define PS4_ITERS0         [ebx + 704]
#define PS4_ITERS1         [ebx + 704 + 4]
#define PS4_ITERS2         [ebx + 704 + 8]
#define PS4_ITERS3         [ebx + 704 + 12]
#define PS4_ITERCTR_L      [ebx + 896]
#define PS4_ITERCTR_H      [ebx + 900]
#define PS4_CUR_MAX_ITERS  [ebx + 920]

// Aliases for 8 single precision points
#define PS8_X03            PS4_X01
#define PS8_X47            PS4_X23
#define PS8_Y03            PS4_Y01
#define PS8_Y47            PS4_Y23
#define PS8_YY03           PS4_YY01
#define PS8_YY47           PS4_YY23
#define PS8_A03            PS4_A01
#define PS8_A47            PS4_A23
#define PS8_B03            PS4_B01
#define PS8_B47            PS4_B23
#define PS8_MAG03          [ebx + 320]
#define PS8_MEXP0          [ebx + 320]
#define PS8_MEXP1          [ebx + 320 + 4]
#define PS8_MEXP2          [ebx + 320 + 8]
#define PS8_MEXP3          [ebx + 320 + 12]
#define PS8_MAG47          [ebx + 320 + 16]
#define PS8_MEXP4          [ebx + 320 + 16]
#define PS8_MEXP5          [ebx + 320 + 20]
#define PS8_MEXP6          [ebx + 320 + 24]
#define PS8_MEXP7          [ebx + 320 + 28]
#define PS8_MAGPREV03      PS4_MAGPREV01
#define PS8_MAGPREV47      PS4_MAGPREV23
#define PS8_TWO            [ebx + 512]
#define PS8_RAD            [ebx + 640]

// Iteration counters for 8 points
#define PS8_ITERS0         [ebx + 704]

#define PS8_ITERS1         [ebx + 704 + 4]
#define PS8_ITERS2         [ebx + 704 + 8]
#define PS8_ITERS3         [ebx + 704 + 12]
#define PS8_ITERS4         [ebx + 704 + 16]
#define PS8_ITERS5         [ebx + 704 + 20]
#define PS8_ITERS6         [ebx + 704 + 24]
#define PS8_ITERS7         [ebx + 704 + 28]
#define PS8_ITERCTR_L      PS4_ITERCTR_L
#define PS8_ITERCTR_H      PS4_ITERCTR_H
#define PS8_CUR_MAX_ITERS  PS4_CUR_MAX_ITERS

/**
 * Generic setting structure. Some flags are encoded by upper or lowercase 
 * letters in the name (see below). The name is not case sensitive in files.
 */
typedef struct {
	char* name; // Name of the setting
	int c_val;	// Value of the setting
	int d_val;	// Default value for the setting
	int min;	// Its min and max limits, to prevent bad data
	int max;	// in files from messing things up
} setting;

#define LOWERCASE 0x20

/**
 * Autoreset flag is set by having the first letter of the name uppercase.
 * If the autoreset flag is set, the setting resets to default_val before each 
 * new image.
 */
#define SETTING_AUTORESET(s) (!((s)->name[0] & LOWERCASE))

/**
 * Global configuration and algorithm settings. 
 * All fields must be setting structs.
 */
typedef struct {
	setting pan_rate;	// value for the GUI pan rate slider
	
	// if nonzero, this image should be auto-panned (contains key code)
	setting pan_key;
	
	setting zoom_rate;	// value for the GUI zoom rate slider
	
	// 0 = no zoom, 1 = zoom in, 2 = zoom out (to be implemented)
	setting zoom_in_out;
	
	// X image size. 0 = maximize window (to be implemented)
	setting xsize;
	
	// Y image size  0 = restore window (to be implemented)
	setting ysize;
	
	// the color for points with max iters- for making lakes and such. Cool...
	setting max_iters_color;
	
	// value to XOR with the palette to create the image palette (0xFFFFFF for 
	// invert, etc)
	setting pal_xor;
	
	setting options;	// options bitfield

	// Stripes per thread (bitfield; 4 bits per num_threads index). For 
	// example, bits 7-4 give the number of stripes per thread for two threads, 
	// 11-8 are for four threads, etc.

	setting stripes;
	
	// thickness of stripes used in striped_blit
	setting blit_stripe_thickness;   
	
	setting pfcmin;		// pan filter constant min and max
	setting pfcmax;		// 10000 times the real value for these
} settings;

/**
 * Default for stripes per thread bitfield.
 * 32 threads: 2; 
 * 16 threads: 3; 
 * 8 threads: 4; 
 * 4 threads: 4; 
 * 2 threads: 7; 
 * 1 thread: 1
 */
#define SPT_DEFAULT  0x234471

/* ---------------------- Bits in options bit-field ------------------------ */

// 1 to recalculate immediately whenever the window is resized, 0 = not
#define OPT_RECALC_ON_RESIZE	1

// 1 if the control dialog should initially be visible after entering 
// fullscreen mode, 0 if not can always toggle on/off with C key.
#define OPT_FULLSCREEN			2

// if 1, start with the normalized rendering algorithm (otherwise standard)
#define OPT_NORMALIZED			4

// if 1, start with an exact algorithm (default fast)
#define OPT_EXACT_ALG			8

// Default for options bit-field
#define OPT_DEFAULT				(OPT_RECALC_ON_RESIZE)

/**
 * A log entry structure. It can have its own set of settings. The settings 
 * fields are all initialized to -1, then set to any values in the logfile that 
 * are found for this entry. Any setting that's >= 0 will be copied to the 
 * global config settings before the image is displayed.
 */
typedef struct {
	double re;
	double im;
	double mag;
	unsigned int max_iters;
	unsigned int palette;

	// Kinda wasteful to have the entire settings struct here when all we need 
	// are the val fields, but it's easier this way
	settings log_settings;
} log_entry;

/**
 * A stripe for the thread function to calculate.
 * See man_calculate().
 */
typedef struct {
	int xstart;
	int xend;
	int ystart;
	int yend;
} stripe;

/**
 * Maximum number of stripes per image to give each thread. 
 * See man_calculate().
 */
#define MAX_STRIPES		8

/**
 * Thread state structure
 */
typedef struct {
	int thread_num; // thread number

	// pointer to this thread's iterating point structure, from array in 
	// man_calc_struct
	man_pointstruct* ps_ptr;
	
	// list of stripes for this thread to calculate
	stripe stripes[MAX_STRIPES];
	
	int num_stripes;
	void* done_event;	// event set by thread when calculation is finished
	void* calc_struct;	// pointer to parent man_calc_struct

	// Nonessential variables (for profiling, load balance testing, etc)
	unsigned long long total_iters;		// iters value that keeps accumulating 
										// until reset (before next zoom, etc)
	unsigned int points_guessed;		// points guessed in fast algorithm
} thread_state;

// Use this for overhead-sensitive functions (but nothing else, as it can
// hog registers). Seems to have a small effect.
#define FASTCALL __fastcall

/**
 * Define this to dump any loaded BMP palette to palarray.c as a compilable 
 * array.
*/
// #define DUMP_USER_PAL

typedef struct {
	// Arrays of dwords: msb is 0, lower 3 bytes are RGB triplets.
	// Can have variable-size palettes, given by size field.
	unsigned int* rgb;
	int size;
} palette;

// Structure defining the work for each palette mapping thread to do

/**
 * Some of these values don't need to be in the structure (can be globals) but 
 * it really doesn't make much difference.
 */
typedef struct {
	void* calc_struct;
	unsigned int* dest;
	unsigned int* src;
	unsigned int xsize;
	unsigned int ysize;
	unsigned int* pal;
	unsigned int pal_size;
	unsigned int max_iters_color;
	int thread_num;
} pal_work;

/**
 * Maximum max_iters for which palette lookup array is used.
 */
#define PAL_LOOKUP_MAX	32768

/**
 * Structure for holding all the parameters, sub-structures, and memory 
 * pointers used for a mandelbrot calculation. Separate structures are used for 
 * the main calculation and for saving images, so both can run in parallel.
 */
typedef struct {
	// This MUST be the first structure member, due to alignment requirements
	man_pointstruct pointstruct_array[MAX_THREADS];

	// These function pointers get set based on the algorithm in use

	// Point queueing function (C/SSE/SSE2/x87)
	void (FASTCALL *queue_point)(void *calc_struct, man_pointstruct *ps_ptr, 
		unsigned int *iters_ptr);

	// Iteration function (C/SSE/SSE2/x87, AMD/Intel). Returns the number of 
	// iterations done per point.
	unsigned int (*mandel_iterate)(man_pointstruct* ps_ptr);

	// Pointer for man_calculate function
	double (*man_calculate)(int xstart, int xend, int ystart, int yend);

	// State structures and events for each thread used in the calculation
	thread_state thread_states[MAX_THREADS];
	HANDLE thread_done_events[MAX_THREADS];

	// Image size and offset parameters
	int xsize;
	int ysize;
	int min_dimension;	// dimension
	int img_size;		// image size

	// To reduce the effects of precision loss at high magnifications, panning 
	// is now tracked as an offset from the re/im point. The former method was 
	// to update re/im on every pan.

	long long pan_xoffs;
	long long pan_yoffs;

	// Calculation parameters
	double re;			// real part of image center coordinate
	double im;			// imaginary part of image center coordinate
	double mag;			// magnification
	unsigned int max_iters;			// maximum iterations to do per point
	unsigned int max_iters_last;	// last max_iters used in a calculation

	// GUI parameters, latched before the calculation starts
	int rendering;				// rendering mode
	int precision;				// precision (user-desired)
	int algorithm;				// algorithm mode
	int threads_i;				// number of threads string index
	int cur_alg;				// current algorithm (can switch during panning)
	int num_threads;			// total number of calculation threads
	int num_stripes;			// stripes number per thread
	int sse_support;			// 1 for SSE, 2 for SSE and SSE2
	int all_recalculated;		// flag indicating a recalculation of 
								// the whole image happened
	int precision_loss;			// 1 if precision loss detected on most recent 
								// calculation
	int screen_xpos;
	int screen_ypos;
	int status;					// general status bitfield, sstat
	double iter_time;			// mandelbrot iteration time only
	double file_tot_time;		// total calculation time since file opened
	quadrant quad[4];			// The 4 quadrant bitmaps (each of size 
								// man_xsize x man_ysize)

	// Dynamically allocated arrays
	double* img_re;		// arrays for holding the RE, IM coordinates
	double* img_im;		// of each pixel in the image

	// for dummy line creation: see alloc_man_mem
	unsigned int* iter_data_start;
	
	// iteration counts for each pixel in the image. Converted to a bitmap by 
	// applying the palette.
	unsigned int* iter_data;
	
	// size of one line of iteration data (needs 2 dummy pixels at the end of 
	// each line)
	int iter_data_line_size;

	float* mag_data;	// magnitude (squared) for each point
	int mag_data_offs;	// byte offset of mag_data from iter_data. Could be 
						// negative; must be int

	unsigned char* png_buffer;		// buffer for data to write to PNG file

	// Palette and rendering related items
	unsigned int palette;			// current palette to use
	unsigned int prev_pal;			// previous palette used
	unsigned int pal_xor;			// for inversion (use 0xFFFFFF)
	unsigned int max_iters_color;	// color of max iters points, from 
									// logfile/cfgfile

	pal_work pal_work_array[MAX_THREADS]; // work for palette mapping threads
	void* pal_events[MAX_THREADS];
	
	// palette lookup table; need one extra entry for max_iters
	unsigned int pal_lookup[PAL_LOOKUP_MAX + 1];

	// Misc items
	unsigned int flags;	// See below

	// Pointer offsets into bitmap, precalculated from above values
	int wave_ptr_offs[7][4];
} man_calc_struct;

/**
 * Values for man_calc_struct flags field
 */
#define FLAG_IS_SAVE          1 // 1 if this is a saving structure, 0 for 
								// normal calculation
#define FLAG_CALC_RE_ARRAY    2 // set to 0 on first row when saving, otherwise 
								// 1- reduces overhead

/**
 * Get the magnitude (squared) corresponding to the iteration count at iter_ptr. 
 * Points to an entry in the mag_data array of a man_calc_struct.
 */
#define MAG(m, iter_ptr) *((float *) ((char *)iter_ptr + m->mag_data_offs))

/* ------------------------------- Prototypes ------------------------------ */

// imagesave.c
int png_save_start(char* file, int width, int height);
int png_save_write_row(unsigned char* row);
int png_save_end(void);

// palettes.c
int init_palettes(double diverged_thresh, 
	man_calc_struct* m, 
	man_calc_struct* s);
int load_palette(FILE* fp, man_calc_struct* m);
int load_palette_from_bmp(FILE* fp, man_calc_struct* m);
int get_palette_rgb_val(int ind, char* str, int length, unsigned int* rgb);
void apply_palette(man_calc_struct* m, unsigned int* dest, unsigned int* src, 
	unsigned int xsz, 
	unsigned int ysz);

// quickman.c
void update_re_im(man_calc_struct* m, long long xoffs, long long yoffs);
void update_re_im_mag(man_calc_struct* m, int zoom_box, int in_outn,
	int x0, int y0, int x1, int y1);
void reset_thread_load_counters(man_calc_struct* m);
void reset_quadrants(man_calc_struct* m);
void man_init(man_calc_struct* m, int cfg_options_val, unsigned int flags);
void man_calculate_quadrants(man_calc_struct* m);
void man_setup(man_calc_struct* m, int xstart, int xend, int ystart, int yend);
void pan_image(man_calc_struct* m, int offs_x, int offs_y);
int fast_wave_alg(man_calc_struct* m, man_pointstruct* ps_ptr, stripe* s);
char* get_image_info(man_calc_struct* m, int update_iters_sec);
double get_re_im_offs(man_calc_struct* m, long long offs);
TIME_UNIT get_timer(void);
double get_seconds_elapsed(TIME_UNIT start_time);
int alloc_man_mem(man_calc_struct* m, int width, int height);
void free_man_mem(man_calc_struct* m);
void set_wave_ptr_offs(man_calc_struct* m);
void update_iters(man_calc_struct* m, int up, int down);
void set_home_image(man_calc_struct* m);

#endif /*_QUICKMAN_H_ */