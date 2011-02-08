/***************************************************************************
 *                         ctatools - SWIG file                            *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 * ----------------------------------------------------------------------- *
 * Usage:                                                                  *
 * swig -c++ -python -Wall ctatools.i                                      *
 ***************************************************************************/
%module ctatools

/* __ Support module _____________________________________________________ */
%include "stl.i"
%include "exception.i"
%include "GApplication.i"

/* __ CTA tools __________________________________________________________ */
%include "ctobssim.i"
%include "ctselect.i"
%include "ctbin.i"
%include "ctlike.i"
