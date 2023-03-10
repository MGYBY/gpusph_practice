/*  Copyright (c) 2021 INGV, EDF, UniCT, JHU

    Istituto Nazionale di Geofisica e Vulcanologia, Sezione di Catania, Italy
    Électricité de France, Paris, France
    Università di Catania, Catania, Italy
    Johns Hopkins University, Baltimore (MD), USA

    This file is part of GPUSPH. Project founders:
        Alexis Hérault, Giuseppe Bilotta, Robert A. Dalrymple,
        Eugenio Rustico, Ciro Del Negro
    For a full list of authors and project partners, consult the logs
    and the project website <https://www.gpusph.org>

    GPUSPH is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    GPUSPH is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with GPUSPH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef RWN04Fr1PERIODIC_H
#define	RWN04Fr1PERIODIC_H

#define PROBLEM_API 1
#include "Problem.h"

class rwN04Fr1Periodic: public Problem {
	private:
		double nd; // normalDepth
        double nv; // normalVelocity
		double wl; // Wavelenghth
		double distAmp; // disturbance amplitude
		double channelSin; // sinTheta
		double grav; //gravity

        double mudRho; // mud density
        double mudMu; // mud consistency
        double powerLawN; // power-law index

        double		lx, ly; // computational domain

		bool		m_usePlanes; // use planes or boundaries
public:
	rwN04Fr1Periodic(GlobalData *);
	void initializeParticles(BufferList &, const uint);
	void fillDeviceMap();
};
#endif	/* STANDINGWAVE2D_H */

