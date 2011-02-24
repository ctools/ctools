/***************************************************************************
 *                    ctselect - CTA data selection tool                   *
 * ----------------------------------------------------------------------- *
 *  copyright (C) 2010-2011 by Jurgen Knodlseder                           *
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
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTSELECT_NAME    "ctselect"
#define CTSELECT_VERSION "00-02-01"


/***********************************************************************//**
 * @class ctselect
 *
 * @brief CTA data selection tool interface defintion.
 ***************************************************************************/
class ctselect : public GApplication  {
public:
    // Constructors and destructors
    ctselect(void);
    explicit ctselect(GObservations obs);
    ctselect(int argc, char *argv[]);
    ctselect(const ctselect& app);
    virtual ~ctselect(void);

    // Operators
    ctselect& operator= (const ctselect& app);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           select_events(GCTAObservation* obs, const std::string& filename);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctselect& app);
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

    // Protected members
    GObservations m_obs;        //!< Observations container
};

#endif /* CTSELECT_HPP */
