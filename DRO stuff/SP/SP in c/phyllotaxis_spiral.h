/*
 * Copyright (c) 2016 by General Electric Company. All Rights Reserved.
 */

/**
 * \file phyllotaxis_spiral.h
 *
 * This is the header file for Cartesian spiral support routines.
 *
 * @author R. Marc Lebel 
 * @since 25.0
 */

/*
 * Comments:
 * Date (dd mmm yyyy) Author (or Initials)
 * Comment
 *
 * 12 May 2016 RML
 *  Initial version.
 */

#ifndef phyllotaxis_spiral_h
#define phyllotaxis_spiral_h

int gen_phyllotaxis_spiral(const char *filename,
                           int yres,
                           int zres,
                           int N,
                           int pts,
                           double target_rot,
                           double density);

#endif
