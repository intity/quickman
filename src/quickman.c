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

#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "quickman.h"

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

/* ------- Constants and variables used in the fast "wave" algorithm ------- */

/**
 * Starting values for x and y. Now seems faster to have these as static 
 * globals calculates y = 0, needs dummy line at y = -1 (set to 0's)
 */
static const int wave_ystart[7] = { 3, 1, 3, 1, 0, 1, 0 };
static const int wave_xstart[7] = { 0, 2, 2, 0, 1, 1, 0 };

/**
 * X and Y increments for each wave
 */
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

/*-------------------------- Benchmark functions ----------------------------*/

/**
 * Reset frames/sec timing values.
 */
void reset_fps_values(man_calc_struct* m)
{
	m->total_frames       = 0;
	m->interval_frames    = 0;
	m->total_time         = 0.0;
	m->interval_time      = 0.0;
	m->calc_total_time    = 0.0;
	m->calc_interval_time = 0.0;
}

/*-------------------------- File/misc functions ----------------------------*/

/**
 * Timer function: returns the current time in real format.
 */
double get_timer(void)
{
	struct timespec ts;
	timespec_get(&ts, TIME_UTC);

	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/**
 * Get the number of seconds elapsed since start_time.
 */
double get_seconds_elapsed(double start_time)
{
	struct timespec ts;
	timespec_get(&ts, TIME_UTC);

	return (ts.tv_sec + ts.tv_nsec * 1e-9) - start_time;
}

/**
 * Reset thread iters accumulators, used to measure thread load balance. Only
 * used for main calculation.
 */
void reset_thread_load_counters(man_calc_struct* m)
{
	int i;

	for (i = 0; i < MAX_THREADS; i++)
		m->thread_states[i].total_iters = 0;
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
void update_re_im_mag(man_calc_struct* m, int zoom_box, int in_outn, 
	int x0, int y0, 
	int x1, int y1)
{
	double tmp_mag, xz, yz;
	int x, y;

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

/**
 * Queue a point to be iterated, for the C iteration function.
 */
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

/* ------------------- Code for ASM iteration functions -------------------- */

/**
 * XMM registers for the SSE2 algorithms
 */
#define xmm_x01		xmm0
#define xmm_y01		xmm1
#define xmm_yy01	xmm2
#define xmm_mag		xmm3
#define xmm_x23		xmm4
#define xmm_y23		xmm5
#define xmm_yy23	xmm6
#define xmm_two		xmm7

/**
 * XMM registers for the SSE algorithms
 */
#define xmm_x03		xmm0
#define xmm_y03		xmm1
#define xmm_yy03	xmm2
#define xmm_x47		xmm4
#define xmm_y47		xmm5
#define xmm_yy47	xmm6

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

		mov		ebx, ps_ptr			// Get pointstruct pointer. PS4_ macros 
									// below reference [ebx + offset]
		movapd	xmm_x23, PS4_X23	// Restore point states
		movapd	xmm_x01, PS4_X01
		movapd	xmm_two, PS4_TWO
		movapd	xmm_y23, PS4_Y23
		movapd	xmm_y01, PS4_Y01
		movapd	xmm_yy23, PS4_YY23

		mov edx, DIV_EXP			// Exp for magnitude exponent comparison. 
									// Slower to compare to const directly
		mov ecx, PS4_CUR_MAX_ITERS	// max iters to do this call; always even
		mov eax, 0			// iteration counter (for each of the 4 points)
		jmp skip_top		// Jump past 1st 2 movapds for loop entry, 
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
		movapd	PS4_YY01, xmm_yy01	// save yy01 for mag backout checking
		movapd	PS4_X01, xmm_x01	// save x01 for mag backout checking; 
									// contains xx01 - yy01 here
	skip_top:
		addpd	xmm_y01, PS4_B01	// y01 += b01; faster at top of loop. 
									// Initial y01 = 0
		mulpd	xmm_y23, xmm_x23	/* y23 *= x23; faster at top of loop. */
		add eax, 2					// update iteration counter; faster here 
									// than 2 insts below
		mulpd	xmm_x23, xmm_x23	/* x23 *= x23 */
		addpd	xmm_x01, PS4_A01	// x01 += a01
		addpd	xmm_y23, xmm_y23	/* y23 *= 2; faster here than mulpd 
									   xmm_y23, xmm_two */
		movapd	xmm_yy01, xmm_y01	// yy01 = y01
		movapd	PS4_X23, xmm_x23	// save xx23 for magnitude backout checking
		mulpd	xmm_yy01, xmm_yy01	/* yy01 *= yy01 */
		subpd	xmm_x23, xmm_yy23	// x23 -= yy23
		mulpd	xmm_y01, xmm_x01	/* y01 *= x01 */
		addpd	xmm_y23, PS4_B23	// y23 += b23
		mulpd	xmm_x01, xmm_x01	/* x01 *= x01 */
		movapd	PS4_YY23, xmm_yy23	// save yy23 for magnitude backout checking
		mulpd	xmm_y01, xmm_two	/* y01 *= 2; add slower here; bb stall */
		addpd	xmm_x23, PS4_A23	// x23 += a23
		movapd	xmm_mag, xmm_x01	// mag01 = x01
		movapd	xmm_yy23, xmm_y23	// yy23 = y23
		subpd	xmm_x01, xmm_yy01	// x01 -= yy01
		mulpd	xmm_yy23, xmm_yy23	/* yy23 *= yy23 */
		addpd	xmm_y01, PS4_B01	// y01 += b01
		mulpd	xmm_y23, xmm_x23	/* y23 *= x23 */
		// ----- Start of 2nd iteration block ------
		addpd	xmm_mag, xmm_yy01
		mulpd	xmm_x23, xmm_x23
		addpd	xmm_x01, PS4_A01
		movapd	xmm_yy01, xmm_y01	// these 2 instrs: faster in this order 
									// than reversed (fixes y23 dep?)
		mulpd	xmm_y23, xmm_two	// (yy01 is apparently just "marked" to get 
									// y01 here; doesn't cause a dep delay)
		movapd	PS4_MAG01, xmm_mag	// save point 0,1 magnitudes for comparison
		movapd	xmm_mag, xmm_x23
		mulpd	xmm_yy01, xmm_yy01
		subpd	xmm_x23, xmm_yy23
		mulpd	xmm_y01, xmm_x01
		addpd	xmm_y23, PS4_B23
		mulpd	xmm_x01, xmm_x01
		addpd	xmm_mag, xmm_yy23
		movapd	PS4_MAG23, xmm_mag	// Save point 2,3 magnitudes. Best here, 
									// despite dep
		cmp		PS4_MEXP0, edx		// Compare the magnitude exponents of 
									// points 0 and 1
		cmovge	ecx, eax			// to the divergence threshold. AMD doesn't 
									// seem to mind
		cmp		PS4_MEXP1, edx	// the store fowarding issue, but Intel does.
		cmovge	ecx, eax			// Conditional moves set ecx to eax on 
									// divergence,
		mulpd	xmm_y01, xmm_two	// breaking the loop.
		addpd	xmm_x23, PS4_A23	// add y, y and mul y, two seem equal speed 
									// here
		movapd	xmm_yy23, xmm_y23
		subpd	xmm_x01, xmm_yy01
		mulpd	xmm_yy23, xmm_yy23
		cmp		PS4_MEXP2, edx		// Compare the magnitude exponents of 
									// points 2 and 3.
		cmovge	ecx, eax
		cmp		PS4_MEXP3, edx
		cmovge	ecx, eax
		cmp		ecx, eax			// Continue iterating until max iters 
									// reached for this call,
		jne		iter_loop			// or one of the points diverged.

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

		movapd	PS4_Y01, xmm_y01	// save y01 state
		movapd	PS4_Y23, xmm_y23	// save y23 state

		mulpd	xmm_two, PS4_YY01	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy01 */
		movapd	xmm_mag, PS4_X23	// tmp2 = xx23
		addpd	xmm_two, PS4_X01	/* get mag01 = xx01 - yy01 + 2 * yy01 = 
									   xx01 + yy01 */
		addpd	xmm_mag, PS4_YY23	// get mag23 = xx23 + yy23
		movapd	PS4_MAGPREV01, xmm_two		// store prev_mag 01
		movapd	PS4_MAGPREV23, xmm_mag		// store prev_mag 23

		xor ecx, ecx				// Get a 0
		add		PS4_ITERCTR_L, eax	// Update iteration counter. Multiply by 36 
									// to get effective flops.
		adc		PS4_ITERCTR_H, ecx	// Update iterctr high dword

		movapd	PS4_YY01, xmm_yy01	// save yy01 state
		movapd	PS4_X23, xmm_x23	// save x23 state
		movapd	PS4_X01, xmm_x01	// save x01 state
		movapd	PS4_YY23, xmm_yy23	// save yy23 state

		add		PS4_ITERS0, eax	// update point iteration counts
		add		PS4_ITERS1, eax
		add		PS4_ITERS2, eax
		add		PS4_ITERS3, eax

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
		mov		ebx, ps_ptr			// Get pointstruct pointer
		movapd	xmm_yy01, PS4_YY01	// Restore point states
		movapd	xmm_y01, PS4_Y01
		movapd	xmm_x23, PS4_X23
		movapd	xmm_x01, PS4_X01
		movapd	xmm_two, PS4_TWO
		movapd	xmm_y23, PS4_Y23
		movapd	xmm_yy23, PS4_YY23
		addpd	xmm_y01, PS4_B01	// pre-add y01 to get correct initial 
									// condition

		mov eax, 2			// iteration counter (for each of the 4 points)
		jmp skip_top

	iter_loop:				// alignment doesn't seem to matter on Intel
		movapd	PS4_YY01, xmm_yy01	// save yy01 for mag backout checking
		movapd	PS4_X01, xmm_x01	// save x01 for mag backout checking; 
									// contains xx01 - yy01 here
		add eax, 2					// update iteration counter
	skip_top:
		mulpd	xmm_x23, xmm_x23	/* x23 *= x23 */
		addpd	xmm_x01, PS4_A01	// x01 += a01
		addpd	xmm_y23, xmm_y23	/* y23 *= 2; faster here than mulpd 
									   xmm_y23, xmm_two */
		movapd	xmm_yy01, xmm_y01	// yy01 = y01
		movapd	PS4_X23, xmm_x23	// save xx23 for magnitude backout checking
		mulpd	xmm_yy01, xmm_yy01	/* yy01 *= yy01 */
		subpd	xmm_x23, xmm_yy23	// x23 -= yy23
		mulpd	xmm_y01, xmm_x01	/* y01 *= x01 */
		addpd	xmm_y23, PS4_B23	// y23 += b23
		mulpd	xmm_x01, xmm_x01	/* x01 *= x01 */
		movapd	PS4_YY23, xmm_yy23	// save yy23 for magnitude backout checking
		mulpd	xmm_y01, xmm_two	/* y01 *= 2; add slower here */
		addpd	xmm_x23, PS4_A23	// x23 += a23
		movapd	xmm_mag, xmm_x01	// mag01 = x01
		movapd	xmm_yy23, xmm_y23	// yy23 = y23
		subpd	xmm_x01, xmm_yy01	// x01 -= yy01
		mulpd	xmm_yy23, xmm_yy23	/* yy23 *= yy23 */
		addpd	xmm_y01, PS4_B01	// y01 += b01
		mulpd	xmm_y23, xmm_x23	/* y23 *= x23 */
		// ----- Start of 2nd iteration block ------
		addpd	xmm_mag, xmm_yy01
		mulpd	xmm_x23, xmm_x23
		addpd	xmm_x01, PS4_A01
		mulpd	xmm_y23, xmm_two
		movapd	xmm_yy01, xmm_y01
		movapd	PS4_MAG01, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltpd xmm_mag, PS4_RAD	// compare point 0, 1 magnitudes 
									// (mag >= rad): let cpu reorder these
		movmskpd edx, xmm_mag		// save result in edx
		movapd	xmm_mag, xmm_x23
		mulpd	xmm_yy01, xmm_yy01
		subpd	xmm_x23, xmm_yy23
		mulpd	xmm_y01, xmm_x01
		addpd	xmm_y23, PS4_B23
		mulpd	xmm_x01, xmm_x01
		addpd	xmm_mag, xmm_yy23
		mulpd	xmm_y01, xmm_two
		addpd	xmm_x23, PS4_A23
		add edx, edx			// shift point 01 mag compare results left 2
		add edx, edx
		movapd	xmm_yy23, xmm_y23
		subpd	xmm_x01, xmm_yy01
		movapd	PS4_MAG23, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltpd xmm_mag, PS4_RAD	// compare point 2, 3 magnitudes
		mulpd	xmm_yy23, xmm_yy23
		addpd	xmm_y01, PS4_B01
		movmskpd ecx, xmm_mag
		mulpd	xmm_y23, xmm_x23
		or ecx,	edx					// Continue iterating until max iters 
									// reached for this call,
		jnz	done					// or one of the points diverged.
		cmp	PS4_CUR_MAX_ITERS, eax	// No penalty for comparing from 
									// memory vs. register here
		jne iter_loop

	done:
		subpd	xmm_y01, PS4_B01	// subtract out pre-add (see loop top)
		movapd	PS4_Y01, xmm_y01	// save y01 state
		movapd	PS4_Y23, xmm_y23	// save y23 state

		// Get previous magnitudes. See AMD code
		mulpd	xmm_two, PS4_YY01	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy01 */
		movapd	xmm_mag, PS4_X23	// tmp2 = xx23
		addpd	xmm_two, PS4_X01	/* get mag01 = xx01 - yy01 + 2 * yy01 = 
									   xx01 + yy01 */
		addpd	xmm_mag, PS4_YY23	// get mag23 = xx23 + yy23
		movapd	PS4_MAGPREV01, xmm_two // store prev_mag 01
		movapd	PS4_MAGPREV23, xmm_mag // store prev_mag 23

		xor edx, edx				// Get a 0
		add		PS4_ITERCTR_L, eax	// Update iteration counter. Multiply by 36 
									// to get effective flops.
		adc		PS4_ITERCTR_H, edx	// Update iterctr high dword

		movapd	PS4_YY01, xmm_yy01	// save yy01 state
		movapd	PS4_X23, xmm_x23	// save x23 state
		movapd	PS4_X01, xmm_x01	// save x01 state
		movapd	PS4_YY23, xmm_yy23	// save yy23 state

		add		PS4_ITERS0, eax		// update point iteration counts
		add		PS4_ITERS1, eax
		add		PS4_ITERS2, eax
		add		PS4_ITERS3, eax
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
		mov		ebx, ps_ptr			// Get pointstruct pointer
		movaps	xmm_x47, PS8_X47	// Restore point states
		movaps	xmm_x03, PS8_X03
		movaps	xmm_two, PS8_TWO
		movaps	xmm_y47, PS8_Y47
		movaps	xmm_y03, PS8_Y03
		movaps	xmm_yy47, PS8_YY47

		mov	edx, DIV_EXP_FLOAT		// Exp for magnitude exponent comparison. 
									// Slower to compare to const directly
		mov ecx, PS8_CUR_MAX_ITERS	// max iters to do this call; always even
		mov eax, 0		// Iteration counter (for each of the 4 points)
		jmp skip_top				// Allows removing yy03 restore above
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
		movaps	PS8_YY03, xmm_yy03	// save yy03 for mag backout checking
		movaps	PS8_X03, xmm_x03	// save x03 for mag backout checking; 
									// contains xx03 - yy03 here
	skip_top:
		addps	xmm_y03, PS8_B03	// y03 += b03; faster at top of loop. 
									// Initial y03 = 0
		mulps	xmm_y47, xmm_x47	/* y47 *= x47; faster at top of loop. */
		add eax, 2					// update iteration counter; faster here 
									// than 2 insts below
		mulps	xmm_x47, xmm_x47	// x47 *= x47
		addps	xmm_x03, PS8_A03	// x03 += a03
		addps	xmm_y47, xmm_y47	/* y47 *= 2; faster here than mulps 
									   xmm_y47, xmm_two */
		movaps	xmm_yy03, xmm_y03	// yy03 = y03
		movaps	PS8_X47, xmm_x47	// save xx47 for magnitude backout checking
		mulps	xmm_yy03, xmm_yy03	/* yy03 *= yy03 */
		subps	xmm_x47, xmm_yy47	// x47 -= yy47
		mulps	xmm_y03, xmm_x03	/* y03 *= x03 */
		addps	xmm_y47, PS8_B47	// y47 += b47
		mulps	xmm_x03, xmm_x03	/* x03 *= x03 */
		movaps	PS8_YY47, xmm_yy47	// save yy47 for magnitude backout checking
		mulps	xmm_y03, xmm_two	/* y03 *= 2; add slower here; bb stall */
		addps	xmm_x47, PS8_A47	// x47 += a47
		movaps	xmm_mag, xmm_x03	// mag03 = x03
		movaps	xmm_yy47, xmm_y47	// yy47 = y47
		subps	xmm_x03, xmm_yy03	// x03 -= yy03
		mulps	xmm_yy47, xmm_yy47	/* yy47 *= yy47 */
		addps	xmm_y03, PS8_B03	// y03 += b03
		mulps	xmm_y47, xmm_x47	/* y47 *= x47 */
		// ----- Start of 2nd iteration block ------
		addps	xmm_mag, xmm_yy03
		mulps	xmm_x47, xmm_x47
		addps	xmm_x03, PS8_A03
		movaps	xmm_yy03, xmm_y03	// these 2 instrs: faster in this order 
									// than reversed (fixes y47 dep?)
		mulps	xmm_y47, xmm_two	// (yy03 is apparently just "marked" to get 
									// y03 here; doesn't cause a dep delay)
		movaps	PS8_MAG03, xmm_mag	// save point 0-3 magnitudes for comparison
		movaps	xmm_mag, xmm_x47
		mulps	xmm_yy03, xmm_yy03
		subps	xmm_x47, xmm_yy47
		mulps	xmm_y03, xmm_x03
		addps	xmm_y47, PS8_B47
		mulps	xmm_x03, xmm_x03
		addps	xmm_mag, xmm_yy47
		movaps	PS8_MAG47, xmm_mag	// Save point 4-7 magnitudes. Best here, 
									// despite dep
		cmp		PS8_MEXP0, edx		// Compare the magnitude exponents of 
									// points 0-3
		cmovge	ecx, eax			// to the divergence threshold. AMD doesn't 
									// seem to mind
		cmp		PS8_MEXP1, edx		// the store fowarding issue, but Intel 
									// does.
		cmovge	ecx, eax			// Conditional moves set ecx to eax on 
									// divergence,
		cmp		PS8_MEXP2, edx		// breaking the loop.
		cmovge	ecx, eax
		cmp		PS8_MEXP3, edx
		cmovge	ecx, eax
		mulps	xmm_y03, xmm_two	// add y, y and mul y, two seem equal speed 
									// here
		addps	xmm_x47, PS8_A47
		movaps	xmm_yy47, xmm_y47
		subps	xmm_x03, xmm_yy03
		mulps	xmm_yy47, xmm_yy47
		cmp		PS8_MEXP4, edx		// Compare the magnitude exponents of 
									// points 4-7.
		cmovge	ecx, eax
		cmp		PS8_MEXP5, edx
		cmovge	ecx, eax
		cmp		PS8_MEXP6, edx
		cmovge	ecx, eax
		cmp		PS8_MEXP7, edx
		//cmovge	ecx, eax			// jge done seems a hair faster. When 
									// changing instrs, don't forget
		jge	done					// to adjust nops to put jne iter_loop on a 
									// 16-byte boundary.
		cmp	ecx, eax				// Continue iterating until max iters 
									// reached for this call,
		jne iter_loop				// or one of the points diverged.

	done:
		// Get previous magnitudes. See AMD SSE2 code
		movaps	PS8_Y03, xmm_y03	// save y03 state
		movaps	PS8_Y47, xmm_y47	// save y47 state

		mulps	xmm_two, PS8_YY03	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy03 */
		movaps	xmm_mag, PS8_X47	// tmp2 = xx47
		addps	xmm_two, PS8_X03	/* get mag03 = xx03 - yy03 + 2 * yy03 = 
									   xx03 + yy03 */
		addps	xmm_mag, PS8_YY47	// get mag47 = xx47 + yy47
		movaps	PS8_MAGPREV03, xmm_two		// store prev_mag 03
		movaps	PS8_MAGPREV47, xmm_mag		// store prev_mag 47

		xor ecx, ecx				// Get a 0
		add		PS8_ITERCTR_L, eax	// Update iteration counter. Multiply by 72 
									// to get effective flops.
		adc		PS8_ITERCTR_H, ecx	// Update iterctr high dword

		movaps	PS8_YY03, xmm_yy03	// save yy03 state
		movaps	PS8_X47, xmm_x47	// save x47 state
		movaps	PS8_X03, xmm_x03	// save x03 state
		movaps	PS8_YY47, xmm_yy47	// save yy47 state

		add		PS8_ITERS0, eax		// update point iteration counts
		add		PS8_ITERS1, eax
		add		PS8_ITERS2, eax
		add		PS8_ITERS3, eax
		add		PS8_ITERS4, eax
		add		PS8_ITERS5, eax
		add		PS8_ITERS6, eax
		add		PS8_ITERS7, eax
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
		mov		ebx, ps_ptr			// Get pointstruct pointer
		movaps	xmm_yy03, PS8_YY03	// Restore point states
		movaps	xmm_x47, PS8_X47
		movaps	xmm_x03, PS8_X03
		movaps	xmm_two, PS8_TWO
		movaps	xmm_y47, PS8_Y47
		movaps	xmm_y03, PS8_Y03
		movaps	xmm_yy47, PS8_YY47
		addps	xmm_y03, PS8_B03	// pre-add y03 to get correct initial 
									//condition
		mov eax, 2			// iteration counter (for each of the 4 points)
		jmp skip_top

	iter_loop:					// alignment doesn't seem to matter on Intel
		movaps	PS8_YY03, xmm_yy03	// save yy03 for mag backout checking
		movaps	PS8_X03, xmm_x03	// save x03 for mag backout checking; 
									// contains xx03 - yy03 here
		add eax, 2					// update iteration counter
	skip_top:
		mulps	xmm_x47, xmm_x47	/* x47 *= x47 */
		addps	xmm_x03, PS8_A03	// x03 += a03
		addps	xmm_y47, xmm_y47	/* y47 *= 2; faster here than mulps 
									   xmm_y47, xmm_two */
		movaps	xmm_yy03, xmm_y03	// yy03 = y03
		movaps	PS8_X47, xmm_x47	// save xx47 for magnitude backout checking
		mulps	xmm_yy03, xmm_yy03	/* yy03 *= yy03 */
		subps	xmm_x47, xmm_yy47	// x47 -= yy47
		mulps	xmm_y03, xmm_x03	/* y03 *= x03 */
		addps	xmm_y47, PS8_B47	// y47 += b47
		mulps	xmm_x03, xmm_x03	/* x03 *= x03 */
		movaps	PS8_YY47, xmm_yy47	// save yy47 for magnitude backout checking
		mulps	xmm_y03, xmm_two	/* y03 *= 2; add slower here */
		addps	xmm_x47, PS8_A47	// x47 += a47
		movaps	xmm_mag, xmm_x03	// mag03 = x03
		movaps	xmm_yy47, xmm_y47	// yy47 = y47
		subps	xmm_x03, xmm_yy03	// x03 -= yy03
		mulps	xmm_yy47, xmm_yy47	/* yy47 *= yy47 */
		addps	xmm_y03, PS8_B03	// y03 += b03
		mulps	xmm_y47, xmm_x47	/* y47 *= x47 */
		// ----- Start of 2nd iteration block ------
		addps	xmm_mag, xmm_yy03
		mulps	xmm_x47, xmm_x47
		addps	xmm_x03, PS8_A03
		mulps	xmm_y47, xmm_two
		movaps	xmm_yy03, xmm_y03
		
		movapd	PS8_MAG03, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltps xmm_mag, PS8_RAD	// compare point 0-3 magnitudes 
									// (mag >= rad): let cpu reorder these
		movmskps edx, xmm_mag		// save result in edx
		movaps	xmm_mag, xmm_x47
		mulps	xmm_yy03, xmm_yy03
		subps	xmm_x47, xmm_yy47
		mulps	xmm_y03, xmm_x03
		addps	xmm_y47, PS8_B47
		mulps	xmm_x03, xmm_x03
		addps	xmm_mag, xmm_yy47
		mulps	xmm_y03, xmm_two
		addps	xmm_x47, PS8_A47
		shl	edx, 4				// shift point 0-3 mag compare results left 4
		movaps	xmm_yy47, xmm_y47
		subps	xmm_x03, xmm_yy03
		movapd	PS8_MAG47, xmm_mag	// new, mag store for normalized iteration 
									// count alg -- not much effect on speed
		cmpnltps xmm_mag, PS8_RAD	// compare point 4-7 magnitudes
		mulps	xmm_yy47, xmm_yy47
		addps	xmm_y03, PS8_B03
		movmskps ecx, xmm_mag
		mulps	xmm_y47, xmm_x47
		or ecx, edx		// Continue iterating until max iters reached for this 
						// call,
		jnz		done	// or one of the points diverged.
		cmp		PS8_CUR_MAX_ITERS, eax		// No penalty for comparing from 
											// memory vs. register here
		jne		iter_loop

	done:
		subps	xmm_y03, PS8_B03	// subtract out pre-add (see loop top)
		movaps	PS8_Y03, xmm_y03	// save y03 state
		movaps	PS8_Y47, xmm_y47	// save y47 state

		// Get previous magnitudes. See AMD SSE2 code
		mulps	xmm_two, PS8_YY03	/* Use xmm_two for tmp var; 
									   tmp1 = 2 * yy03 */
		movaps	xmm_mag, PS8_X47	// tmp2 = xx47
		addps	xmm_two, PS8_X03	/* get mag03 = xx03 - yy03 + 2 * yy03 = 
									   xx03 + yy03 */
		addps	xmm_mag, PS8_YY47	// get mag47 = xx47 + yy47
		movaps	PS8_MAGPREV03, xmm_two	// store prev_mag 03
		movaps	PS8_MAGPREV47, xmm_mag	// store prev_mag 47

		xor edx, edx				// Get a 0
		add		PS8_ITERCTR_L, eax	// Update iteration counter. Multiply 
									// by 72 to get effective flops.
		adc		PS8_ITERCTR_H, edx	// Update iterctr high dword

		movaps	PS8_YY03, xmm_yy03	// save yy03 state
		movaps	PS8_X47, xmm_x47	// save x47 state
		movaps	PS8_X03, xmm_x03	// save x03 state
		movaps	PS8_YY47, xmm_yy47	// save yy47 state

		add		PS8_ITERS0, eax		// update point iteration counts
		add		PS8_ITERS1, eax
		add		PS8_ITERS2, eax
		add		PS8_ITERS3, eax
		add		PS8_ITERS4, eax
		add		PS8_ITERS5, eax
		add		PS8_ITERS6, eax
		add		PS8_ITERS7, eax
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
#define QUEUE_FULL	0xF

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
	unsigned int* iters_ptr)
{
	unsigned int i, iters, max, queue_status, * ptr;
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
 * Fast "wave" algorithm from old code: guesses pixels.
 */
int fast_wave_alg(man_calc_struct* m, man_pointstruct* ps_ptr, stripe* s)
{
	unsigned int* iters_ptr;
	int x, y, xoffs, offs0, offs1, offs2, offs3;
	int line_size, points_guessed, wave, inc, p0, p1, p2, p3;

	line_size = m->iter_data_line_size;
	points_guessed = 0;

	// Doing the full calculation (all waves) on horizontal chunks to improve 
	// cache locality gives no speedup (tested before realtime zooming was 
	// implemented- maybe should test again).

	for (wave = 0; wave < 7; wave++)
	{
		inc = wave_inc[wave];
		y = wave_ystart[wave] + s->ystart;

		// Special case for wave 0 (always calculates all pixels). 
		// Makes realtime zooming measurably faster. X starts at xstart for 
		// wave 0, so can use do-while. For Y, need to calculate all waves even 
		// if out of range, because subsequent waves look forward to pixels 
		// calculated in previous waves (wave 0 starts at y = 3).

		// it's faster with the special case inside the wave loop than 
		// outside
		if (!wave)
		{
			do
			{
				x = s->xstart;

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
				} while (x <= s->xend);
			} while ((y += inc) <= s->yend);
		}
		else  // waves 1-6 check neighboring pixels
		{
			// pointer offsets of neighboring pixels
			offs0 = m->wave_ptr_offs[wave][0];
			offs1 = m->wave_ptr_offs[wave][1];
			offs2 = m->wave_ptr_offs[wave][2];
			offs3 = m->wave_ptr_offs[wave][3];

			xoffs = wave_xstart[wave] + s->xstart;

			do
			{
				x = xoffs;
				ps_ptr->ab_in[1] = m->img_im[y];
				iters_ptr = m->iter_data + y * line_size + x;

				// No faster to have a special case for waves 1 and 4 
				// that loads only 2 pixels/loop
				while (x <= s->xend)
				{
					// If all 4 neighboring pixels (p0 - p3) are the 
					// same, set this pixel to their value, else 
					// iterate.

					p0 = iters_ptr[offs0];
					p1 = iters_ptr[offs1];
					p2 = iters_ptr[offs2];
					p3 = iters_ptr[offs3];

					// can't use sum compares here (causes corrupted pixels)
					if (p0 == p1 && p0 == p2 && p0 == p3)
					{
						// aargh... compiler (or AMD CPU) generates different 
						// performance on zoomtest depending on which point is 
						// stored here. They're all the same...
						// p3: 18.5s  p2: 18.3s  p1: 18.7s  p0: 19.2s 
						// (+/- 0.1s repeatability)

						*iters_ptr = p2;

						// This works suprisingly well- degradation is really 
						// only noticeable at high frequency transitions (e.g. 
						// with striped palettes).
						// Maybe average the mags at the 4 offsets to make it 
						// better.

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
			} while ((y += inc) <= s->yend);
		}
		// really should flush at the end of each wave, but any errors should 
		// have no visual effect
	}

	return points_guessed;
}

/**
 * Check for precision loss- occurs if the two doubles (or converted floats)
 * in ptest are equal to each other. Returns PLOSS_FLOAT for float loss,
 * PLOSS_DOUBLE for double loss, etc. or 0 for no loss.
 *
 * New more conservative version demands that bits beyond the lsb should also
 * differ. If only the lsb differs, bound to get degradation during iteration.
 */

#define PLOSS_FLOAT		1
#define PLOSS_DOUBLE	2

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
	unsigned int queue_init;
	man_pointstruct* ps_ptr;

	m->max_iters &= ~1;	// make max iters even- required by optimized alg
	xsize = m->xsize;
	ysize = m->ysize;

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

	if (!(m->flags & FLAG_IS_SAVE))
	{
		xend += 4;  // only need these for non-save (precision loss checking)
		yend += 4;
	}

	// this flag should be 1 for main calculation
	if (m->flags & FLAG_CALC_RE_ARRAY)
	{
		ploss = 0;
		x = xstart;
		step = -(xsize >> 1) + xstart + m->pan_xoffs;
		do
		{
			m->img_re[x] = m->re + get_re_im_offs(m, step++);
			if (!(m->flags & FLAG_IS_SAVE) && x > xstart)
				ploss |= check_precision_loss(&m->img_re[x - 1]);
		} while (++x <= xend);
	}

	step = -(ysize >> 1) + ystart + m->pan_yoffs;
	y = ystart;
	do
	{
		m->img_im[y] = m->im - get_re_im_offs(m, step++);
		if (!(m->flags & FLAG_IS_SAVE) && y > ystart)
			ploss |= check_precision_loss(&m->img_im[y - 1]);
	} while (++y <= yend);

	if (!(m->flags & FLAG_IS_SAVE)) // only do auto precision if not saving
	{
		m->precision_loss = 0;

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
				m->precision_loss = 1;
			break;
		case PRECISION_SINGLE:
			if (ploss & PLOSS_FLOAT)
				m->precision_loss = 1;
			break;
		default:
			// should never get here (x87 is suppressed until implemented)
			break;
		}
	}

	// Set iteration and queue_point function pointers and initialize queues

	// Alg will always be C if no sse support.
	// Should change algorithm in dialog box if it's reset to C here.
	if ((m->algorithm & ALG_C) ||
		(m->sse_support < 2 && m->precision == PRECISION_DOUBLE))
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
			m->mandel_iterate = (m->algorithm & ALG_INTEL)
				? iterate_intel_sse2 : iterate_amd_sse2;
		}
		else
		{
			queue_init = (QUEUE_FULL << 24) |
				(7 << 21) | (6 << 18) | (5 << 15) |
				(4 << 12) | (3 << 9) | (2 << 6) | (1 << 3) | 0;

			m->queue_point = queue_8point_sse;
			m->mandel_iterate = (m->algorithm & ALG_INTEL)
				? iterate_intel_sse : iterate_amd_sse;
		}
	}

	if (m->algorithm & ALG_FAST_C)
	{
		queue_init = 0; // must define a default value
	}

	// Set pointstruct initial values
	for (i = 0; i < m->num_threads; i++)
	{
		ps_ptr = m->thread_states[i].ps_ptr;
		ps_ptr->queue_status = queue_init;
		ps_ptr->cur_max_iters = m->max_iters;
		ps_ptr->iterctr = 0;
	}
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
	void* tmp_handle;

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
void reset_quadrants(man_calc_struct* m)
{
	int xsize, ysize;
	quadrant* q;

	xsize = m->xsize;
	ysize = m->ysize;

	q = &m->quad[UL];			// upper left
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

	q = &m->quad[UR];			// upper right
	q->status = 0;				// inactive
	q->quad_rect.x[0] = xsize;	// rectangle coordinates
	q->quad_rect.y[0] = 0;
	q->quad_rect.x[1] = (xsize << 1) - 1;
	q->quad_rect.y[1] = ysize - 1;

	q = &m->quad[LL];			// lower left
	q->status = 0;				// inactive
	q->quad_rect.x[0] = 0;		// rectangle coordinates
	q->quad_rect.y[0] = ysize;
	q->quad_rect.x[1] = xsize - 1;
	q->quad_rect.y[1] = (ysize << 1) - 1;

	q = &m->quad[LR];			// lower right
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

	m->screen_xpos = 0;		// reset screen window to cover UL quadrant only
	m->screen_ypos = 0;
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
void man_calculate_quadrants(man_calc_struct* m)
{
	rectangle r;
	int i, j, x, y;
	unsigned int* bmp_ptr, * iters_ptr;

	m->iter_time = 0.0;

	// First calculate the update rectangles (up to 2).
	for (i = 0; i < 2; i++)
		if (update_rect[i].valid)
		{
			// To get position in (screen-mapped) image, subtract screen upper 
			// left coordinates, Rectangles will be at one of the screen edges 
			// (left, right, top, or bottom).
			// Could simplify this: determined solely by pan offs_x and offs_y

			// Iterate on the update rectangles
			m->iter_time += m->man_calculate(
				update_rect[i].x[0] - m->screen_xpos,	// xstart
				update_rect[i].x[1] - m->screen_xpos,	// xend
				update_rect[i].y[0] - m->screen_ypos,	// ystart
				update_rect[i].y[1] - m->screen_ypos);	// yend
		}

	// Now palette-map the update rectangles into their quadrants. Each 
	// rectangle can occupy 1-4 quadrants.

	for (i = 0; i < 4; i++)
		for (j = 0; j < 2; j++)
			if (update_rect[j].valid)
				if (intersect_rect(&r, &m->quad[i].quad_rect, &update_rect[j]))
				{
					// Subtract upper left coordinates of this quadrant from 
					// upper left coords of intersection rect to get the x, y 
					// offset of the bitmap data in this quadrant.

					x = r.x[0] - m->quad[i].quad_rect.x[0];
					y = r.y[0] - m->quad[i].quad_rect.y[0];

					// get pointer to bitmap data
					bmp_ptr = m->quad[i].bitmap_data + y * m->xsize + x;

					// Subtract upper left coords of screen pos from upper left 
					// coords of intersection rect to get iter data offset

					x = r.x[0] - m->screen_xpos;
					y = r.y[0] - m->screen_ypos;

					// get pointer to iter data to be mapped
					iters_ptr = m->iter_data + y * m->iter_data_line_size + x;

					// Xsize, ysize = rectangle edge lengths
					m->apply_palette(bmp_ptr, iters_ptr, 
						r.x[1] - r.x[0] + 1, 
						r.y[1] - r.y[0] + 1);
				}
}

/**
 * Pan the image by offs_x and offs_y. Sets iter_time (from iteration functions)
 * to the time to do the iteration only.
 */
void pan_image(man_calc_struct* m, int offs_x, int offs_y)
{
	int xsize, ysize, swap_x, swap_y, tmp;
	quadrant* q;
	rectangle* u;

	if (offs_x | offs_y) // Recalculate only if the image moved
	{
		// Update pan offsets
		m->pan_xoffs -= offs_x;	// maybe invert offs_x, offs_y signs later
		m->pan_yoffs -= offs_y;

		xsize = m->xsize;
		ysize = m->ysize;

		// See algorithm explanation above
		m->screen_xpos -= offs_x;  // update screen pos
		m->screen_ypos -= offs_y;

		// Renormalize screen coordinates and swap quadrants if necessary
		swap_x = swap_y = 0;

		if (m->screen_xpos < 0)
		{
			m->screen_xpos += xsize;
			swap_x = 1;
		}
		if (m->screen_xpos > xsize)
		{
			m->screen_xpos -= xsize;
			swap_x = 1;
		}
		if (m->screen_ypos < 0)
		{
			m->screen_ypos += ysize;
			swap_y = 1;
		}
		if (m->screen_ypos > ysize)
		{
			m->screen_ypos -= ysize;
			swap_y = 1;
		}
		if (swap_x)
		{
			swap_quadrants(&m->quad[UL], &m->quad[UR]);
			swap_quadrants(&m->quad[LL], &m->quad[LR]);
		}
		if (swap_y)
		{
			swap_quadrants(&m->quad[UL], &m->quad[LL]);
			swap_quadrants(&m->quad[UR], &m->quad[LR]);
		}

		/* Get the update rectangles (there can be either 1 or 2): determined 
		   by difference between screen pos and previous screen pos (offs_x and 
		   offs_y). Use quadrant coordinates (0,0 to 2 * xsize -1, 2 * ysize -1). 
		   Rectangles will be split at quadrant boundaries during palette 
		   mapping. If either of the sizes is 0, the rectangle will be ignored 
		   (not calculated or palette mapped). */

		// Check for vertical rectangles (have an x offset)
		u = &update_rect[1];
		u->valid = offs_x;						// valid if any nonzero offset
		u->y[0] = m->screen_ypos;				// vertical rectangles are
		u->y[1] = m->screen_ypos + ysize - 1;	// full height of screen

		if (offs_x > 0)
		{  // vertical rect at left
			u->x[0] = m->screen_xpos;
			u->x[1] = m->screen_xpos + offs_x - 1;
		}
		if (offs_x < 0)
		{  // vertical rect at right
			u->x[0] = m->screen_xpos + xsize + offs_x;
			u->x[1] = m->screen_xpos + xsize - 1;
		}

		// Check for horizontal rectangles (have a y offset).
		// Optimize out the corner intersection with any vertical rectangles so 
		// we don't calculate it twice. Clip it off the vertical rect. The 
		// intersection could be a large area for drag pans.

		u = &update_rect[0];
		u->valid = offs_y;						// valid if any nonzero offset
		u->x[0] = m->screen_xpos;				// horizontal rectangles are
		u->x[1] = m->screen_xpos + xsize - 1;	// full width of screen

		if (offs_y > 0)
		{
			// horizontal rect at top
			u->y[0] = m->screen_ypos;
			u->y[1] = tmp = m->screen_ypos + offs_y - 1;
			
			// clip off corner intersection from any vertical rect
			u[1].y[0] = tmp + 1;
		}

		if (offs_y < 0)
		{
			// horizontal rect at bottom
			u->y[0] = tmp = m->screen_ypos + ysize + offs_y;
			u->y[1] = m->screen_ypos + ysize - 1;
			
			// clip off corner intersection from any vertical rect
			u[1].y[1] = tmp - 1;
		}

		// Get the blit rectangles from the screen position (screen_xpos, 
		// screen_ypos = screen upper left corner). Screen coordinates always 
		// range from 0 to xsize and 0 to ysize inclusive here.
		// Refer to diagram.

		m->quad[UL].status = 0; // Default: all inactive
		m->quad[UR].status = 0;
		m->quad[LL].status = 0;
		m->quad[LR].status = 0;

		// Check if UL has a blit rectangle
		if (m->screen_xpos < xsize && m->screen_ypos < ysize)
		{
			q = &m->quad[UL];
			q->status = QSTAT_DO_BLIT; // need a blit
			q->dest_xoffs = 0; // always blits to screen upper left corner
			q->dest_yoffs = 0;
			q->src_xoffs = m->screen_xpos;
			q->src_yoffs = m->screen_ypos;
			q->blit_xsize = xsize - m->screen_xpos;
			q->blit_ysize = ysize - m->screen_ypos;
		}
		// Check if UR has a blit rectangle
		if (m->screen_xpos > 0 && m->screen_ypos < ysize)
		{
			q = &m->quad[UR];
			q->status = QSTAT_DO_BLIT; // need a blit
			q->dest_xoffs = xsize - m->screen_xpos;
			q->dest_yoffs = 0;	// always blits to screen upper edge
			q->src_xoffs = 0;	// always blits from bitmap left edge
			q->src_yoffs = m->screen_ypos;
			q->blit_xsize = m->screen_xpos;
			q->blit_ysize = ysize - m->screen_ypos;
		}
		// Check if LL has a blit rectangle
		if (m->screen_xpos < xsize && m->screen_ypos > 0)
		{
			q = &m->quad[LL];
			q->status = QSTAT_DO_BLIT;	// need a blit
			q->dest_xoffs = 0;			// always blits to screen left edge
			q->dest_yoffs = ysize - m->screen_ypos;
			q->src_xoffs = m->screen_xpos;
			q->src_yoffs = 0;			// always blits from bitmap top edge
			q->blit_xsize = xsize - m->screen_xpos;
			q->blit_ysize = m->screen_ypos;
		}
		// Check if LR has a blit rectangle
		if (m->screen_xpos > 0 && m->screen_ypos > 0)
		{
			q = &m->quad[LR];
			q->status = QSTAT_DO_BLIT;	// need a blit
			q->dest_xoffs = xsize - m->screen_xpos;
			q->dest_yoffs = ysize - m->screen_ypos;
			q->src_xoffs = 0;	// always blits from bitmap upper left corner
			q->src_yoffs = 0;
			q->blit_xsize = m->screen_xpos;
			q->blit_ysize = m->screen_ypos;
		}

		m->status |= STAT_RECALC_FOR_PALETTE;
	}
}

/**
 * Updates image information.
 * If update_iters_sec is 0, won't update iters/sec and gflops (use during 
 * panning, when calculations will be inaccurate due to short calculation times)
 */
void get_image_info(man_calc_struct* m, int update_iters_sec)
{
	int i, points_guessed;
	unsigned long long ictr_raw;
	unsigned long long ictr_total_raw;
	thread_state* t;

	ictr_raw = 0;
	ictr_total_raw = 0;
	points_guessed = 0;

	for (i = 0; i < m->num_threads; i++)
	{
		t = &m->thread_states[i];
		ictr_raw       += t->ps_ptr->iterctr;	// N iterations per tick
		ictr_total_raw += t->total_iters;
		points_guessed += t->points_guessed;
	}

	if (update_iters_sec)
	{
		// For the C versions, each tick is 1 iteration.
		// For the SSE2 ASM versions, each tick is 4 iterations.
		// For the SSE ASM versions, each tick is 8 iterations.

		m->ictr = ictr_raw;

		if (!(m->algorithm & ALG_C))
		{
			if (m->precision == PRECISION_DOUBLE && m->sse_support >= 2)
				m->ictr <<= 2;
			if (m->precision == PRECISION_SINGLE && m->sse_support >= 1)
				m->ictr <<= 3;
		}

		// Prevent division by 0. If the time is in this neighborhood the 
		// iters/sec won't be accurate anyway.
		if (m->iter_time < 0.001)
			m->iter_time = 0.001;

		m->ictr_m      = (double)m->ictr * 1e-6 / m->iter_time;
		m->ictr_g      = m->ictr_m * 9.0 * 1e-3;
		m->avg_iters   = (double)m->ictr / (double)m->img_size;
		m->guessed_pct = 100.0 * (double)points_guessed / (double)m->img_size;

		// Since one flop is optimized out per 18 flops in the ASM versions,
		// factor should really be 8.5 for those. But actually does 9 
		// "effective" flops per iteration.
	}

	// Get each thread's percentage of the total load, to check balance.

	m->max_cur_pct = 0.0;
	m->max_tot_pct = 0.0;

	for (i = 0; i < m->num_threads; i++)
	{
		t = &m->thread_states[i];
		t->cur_pct = (double)t->ps_ptr->iterctr / (double)ictr_raw * 100.0;
		t->tot_pct = (double)t->total_iters / (double)ictr_total_raw * 100.0;

		if (t->cur_pct > m->max_cur_pct)
			m->max_cur_pct = t->cur_pct;

		if (t->tot_pct > m->max_tot_pct)
			m->max_tot_pct = t->tot_pct;
	}

	// Figure of merit: percentage of best possible speed, which occurs with
	// perfect thread load balancing.

	m->max_cur_pct = 100.0 * 100.0 / (m->num_threads * m->max_cur_pct);
	m->max_tot_pct = 100.0 * 100.0 / (m->num_threads * m->max_tot_pct);
}

/* -------------------------- GUI/misc functions --------------------------- */

/**
 * Allocate all the memory needed by the calculation engine. This needs to be 
 * called (after freeing the previous mem) whenever the image size changes.
 */
int alloc_man_mem(man_calc_struct* m, int width, int height)
{
	int n;

	m->iter_data_line_size = width + 2;
	m->img_size = width * height; // new image size

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

/**
 * Free all memory allocated above.
 */
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
 * Precalculate pointer offsets of neighboring pixels for the fast "wave" 
 * algorithm (these only change when image width changes). Not needed for save.
 */
void set_wave_ptr_offs(man_calc_struct* m)
{
	int i, j, sz = m->iter_data_line_size;

	for (i = 1; i < 7; i++)
		for (j = 0; j < 4; j++)
			m->wave_ptr_offs[i][j] = wave_yoffs[i][j] * sz + wave_xoffs[i][j];
}

/**
 * Increase, decrease, or just clip the max iterations.
 */
void update_iters(man_calc_struct* m, int up, int down)
{
	if (up)
		m->max_iters <<= 1;
	if (down)
		m->max_iters >>= 1;
	if (m->max_iters < MIN_ITERS)
		m->max_iters = MIN_ITERS;
	if (m->max_iters > MAX_ITERS)
		m->max_iters = MAX_ITERS;
}

/**
 * Reset mandelbrot parameters to the home image
 */
void set_home_image(man_calc_struct* m)
{
	m->pan_xoffs	= 0;	// Reset any pan offsets
	m->pan_yoffs	= 0;
	m->re			= HOME_RE;
	m->im			= HOME_IM;
	m->mag			= HOME_MAG;
	m->max_iters	= HOME_MAX_ITERS;	// Better to reset the max iters here. 
										// Don't want large #
	update_iters(m, 0, 0);	// from previous image
}
