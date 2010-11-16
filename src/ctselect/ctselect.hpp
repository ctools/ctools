/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010 by Jurgen Knodlseder                                *
 * ----------------------------------------------------------------------- *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
/**
 * @file ctselect.hpp
 * @brief CTA data selection tool definition
 * @author J. Knodlseder
 */

#ifndef CTSELECT_HPP
#define CTSELECT_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"

/* __Definitions _________________________________________________________ */
#define CTSELECT_NAME    "ctselect"
#define CTSELECT_VERSION "v1r0p0"


/***********************************************************************//**
 * @class ctselect
 *
 * @brief CTA data selection tool interface defintion.
 ***************************************************************************/
class ctselect : public GApplication  {
public:
    // Constructors and destructors
    ctselect(void);
    ctselect(int argc, char *argv[]);
    ~ctselect(void);

    // Methods
    void run(void);
    void get_parameters(void);
    void select(void);
    void append_gti(void);
    void copy(void);

protected:
    // Protected methods
    void init_members(void);
    void free_members(void);

    // User parameters
    std::string   m_infile;     //!< Input event list
    std::string   m_outfile;    //!< Output event list
    double        m_ra;         //!< RA of ROI centre
    double        m_dec;        //!< DEC of ROI centre
    double        m_rad;        //!< ROI radius
    double        m_tmin;       //!< Start time
    double        m_tmax;       //!< Stop time
    double        m_emin;       //!< Lower energy
    double        m_emax;       //!< Upper energy
};

#endif /* CTSELECT_HPP */
