/***************************************************************************
 *                      ctbin - CTA data binning tool                      *
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
 * @file ctbin.hpp
 * @brief CTA data binning tool definition
 * @author J. Knodlseder
 */

#ifndef CTBIN_HPP
#define CTBIN_HPP

/* __ Includes ___________________________________________________________ */
#include "GammaLib.hpp"
#include "GCTALib.hpp"

/* __Definitions _________________________________________________________ */
#define CTBIN_NAME    "ctbin"
#define CTBIN_VERSION "00-02-01"


/***********************************************************************//**
 * @class ctbin
 *
 * @brief CTA data binning tool interface defintion
 *
 * This class bins CTA events into a counts map.
 ***************************************************************************/
class ctbin : public GApplication  {
public:
    // Constructors and destructors
    ctbin(void);
    explicit ctbin(GObservations obs);
    ctbin(int argc, char *argv[]);
    ctbin(const ctbin& app);
    virtual ~ctbin(void);

    // Operators
    ctbin& operator= (const ctbin& app);

    // Methods
    void           clear(void);
    void           execute(void);
    void           run(void);
    void           save(void);
    GObservations& obs(void) { return m_obs; }
    void           get_parameters(void);
    void           bin_events(GCTAObservation* obs);

protected:
    // Protected methods
    void init_members(void);
    void copy_members(const ctbin& app);
    void free_members(void);

    // User parameters
    std::string   m_evfile;     //!< Input event list
    std::string   m_outfile;    //!< Output counts map
    double        m_emin;       //!< Lower energy
    double        m_emax;       //!< Upper energy
    int           m_enumbins;   //!< Number of energy bins
    std::string   m_proj;       //!< WCS projection
    std::string   m_coordsys;   //!< Coordinate system
    double        m_xref;       //!< Longitude reference coordinate
    double        m_yref;       //!< Latitude reference coordinate
    double        m_binsz;      //!< Pixel size
    int           m_nxpix;      //!< Number of pixels in longitude
    int           m_nypix;      //!< Number of pixels in latitude

    // Protected members
    GObservations m_obs;        //!< Observation container
};

#endif /* CTBIN_HPP */
