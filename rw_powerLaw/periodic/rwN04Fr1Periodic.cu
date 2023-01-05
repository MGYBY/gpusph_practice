/*  Copyright (c) 2011-2019 INGV, EDF, UniCT, JHU

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

#include <iostream>

#include "rwN04Fr1Periodic.h"
#include "GlobalData.h"
#include "cudasimframework.cu"

rwN04Fr1Periodic::rwN04Fr1Periodic(GlobalData *_gdata) : Problem(_gdata)
{
	// *** user parameters from command line
	// const DensityDiffusionType RHODIFF = get_option("density-diffusion", DELTA_SPH);
	const DensityDiffusionType RHODIFF = get_option("density-diffusion", BREZZI);
	// density diffusion terms: 0 none, 1 Ferrari, 2 Molteni & Colagrossi, 3 Brezzi
	// const int RHODIFF = get_option("density-diffusion", 3);
	const uint ppH = get_option("ppH", 25);
	const bool USE_CCSPH = get_option("use_ccsph", true);

	// use planes in general
	const bool use_planes = get_option("use_planes", false);
	// use a plane for the bottom
	const bool use_bottom_plane = get_option("bottom-plane", use_planes);
	if (use_bottom_plane && !use_planes)
		throw std::invalid_argument("cannot use bottom plane if not using planes");

	// *** Framework setup
	SETUP_FRAMEWORK(
		space_dimensions<R2>,
		periodicity<PERIODIC_X>,
		turbulence_model<LAMINAR_FLOW>,
		// computational_visc<DYNAMIC>,
		// viscosity<DYNAMICVISC>,
		rheology<POWER_LAW>,
		// boundary<DUMMY_BOUNDARY>
		boundary<LJ_BOUNDARY>,
		add_flags<ENABLE_REPACKING>
		// add_flags<ENABLE_INTERNAL_ENERGY>
	).select_options(
		RHODIFF,
		use_planes, add_flags<ENABLE_PLANES>(),
		USE_CCSPH, add_flags<ENABLE_CCSPH>()
	);

	const int mlsIters = get_option("mls",
		(simparams()->densitydiffusiontype != DENSITY_DIFFUSION_NONE) ? 0 : 10);
	if (mlsIters > 0)
		addFilter(MLS_FILTER, mlsIters);

	addPostProcess(SURFACE_DETECTION);

	nd = 0.00209955189; // normal depth
	nv = 0.143385832; // normal velocity
	wl = nd*48.911; // wavelength

	distAmp = 0.225;

	lx = wl;
	ly = 3.0*nd;

	grav = 9.81;

	channelSin = 0.06;
	mudRho = 1120.0;
	powerLawN = 0.40;
	mudMu = 0.14;

	set_deltap(nd/ppH);

	// m_size = make_double3(lx, ly, 1.0);

	setMaxFall(ly);

	// SPH parameters
	simparams()->dtadaptfactor = 0.1;
	simparams()->tend = 32.0;
	simparams()->buildneibsfreq = 10;
	// simparams()->ferrariLengthScale = H;
	simparams()->densityDiffCoeff = 0.05f;

	// Repacking options
	simparams()->repack_maxiter = 10;
	simparams()->repack_a = 0.1;
	simparams()->repack_alpha = 0.1;

	// dyn_thickness = 3.0*m_deltap;

	// Physical parameters
	set_gravity(make_float3(channelSin*grav, (-1.0)*(pow((1-channelSin*channelSin),0.50))*grav, 0.0));
	// g = get_gravity_magnitude();
	// purely for cosmetic reason, let's round the soundspeed to the next
	// integer
	const float c0 = 100.0*nv;
	// auto mud = add_fluid(mudRho);
	add_fluid(mudRho);
	// Surge speed
	setMaxParticleSpeed(5.0*nv);
	// set_equation_of_state(0, 7.0f, c0);
	set_equation_of_state(0, 7.0f, NAN);

	set_consistency_index(0, (mudMu*pow(pow(2.0,0.5),(powerLawN-1.0))));
	set_visc_power_law(0, powerLawN);
	// physparams()->artvisccoeff = 1e-6*10.0/(physparams()->sscoeff[0]*simparams()->slength);

	// Drawing and saving times
	add_writer(VTKWRITER, 0.01);
	add_writer(COMMONWRITER, 0.01);

	// Name of problem used for directory creation
	m_name = "rwN04Fr1PeriodicMini";

	setFillingMethod(Object::BORDER_TANGENT);

	// Building the geometry
	setPositioning(PP_CORNER);
	
	// const int num_layers = (simparams()->boundarytype > SA_BOUNDARY) ?
	// 	simparams()->get_influence_layers() : 1;
	// const double box_thickness = (num_layers - 1)*m_deltap;
	// setDynamicBoundariesLayers(num_layers);
	// // place the walls: as planes, if required; otherwise, as boxes
	// if (use_planes) {
	// 	addPlane(0.0, 1.0, 0.0, 0.0); //bottom plane
	// } else {
	// 	// flat bottom rectangle (before the slope begins)
	// 	GeometryID bottom = addRect(GT_FIXED_BOUNDARY, FT_SOLID,
	// 		Point(paddle_origin - make_double3(box_thickness, box_thickness, 0)),
	// 		h_length + box_thickness + rot_correction, box_thickness);
	// 	setUnfillRadius(bottom, 0.5*m_deltap);
	// }

	const double half_dp = 0.5*m_deltap;
	GeometryID domain_box = addRect(GT_FIXED_BOUNDARY, FT_SOLID,
			Point(half_dp, -3.5*m_deltap, 0), lx-m_deltap, 3*m_deltap);

	// const double half_dp = 0.5*m_deltap;
	// GeometryID domain_box = addRect(GT_FIXED_BOUNDARY, FT_SOLID,
			// Point(half_dp, -3.5*m_deltap, 0), l-m_deltap, 3*m_deltap);

	double3 m_fluidOrigin = make_double3(half_dp, half_dp, 0.0);

	GeometryID fluid = addRect(GT_FLUID, FT_SOLID,
		m_fluidOrigin, lx - m_deltap, nd);
	// double3 m_fluidOrigin = make_double3(0.0, 0.0, 0.0);

	// for now. don't know how to create a complicated free-surface IC
	// GeometryID fluidBox = addRect(GT_FLUID, FT_SOLID,
	// 	m_fluidOrigin, lx, nd);
}

void rwN04Fr1Periodic::initializeParticles(BufferList &buffer, const uint numParticle)
	{


		double4 *gpos = buffer.getData<BUFFER_POS_GLOBAL>();
		float4 *pos = buffer.getData<BUFFER_POS>();
		float4 *vel = buffer.getData<BUFFER_VEL>();
		const ushort4 *pinfo = buffer.getData<BUFFER_INFO>();

		for (uint i = 0 ; i < numParticle ; i++) {
			if (FLUID(pinfo[i])){

				double4 pg = gpos[i];
				// vel[i].x = nv*(1.0+distAmp*sin(2.0*M_PI*pg.x/wl));
				// vel[i].x = (1.0+2.0*powerLawN)/(1.0+powerLawN)*(nv*(1.0+distAmp*sin(2.0*M_PI*pg.x/wl)))*(1.0-(pow((1.0-pg.y/nd),((1.0+powerLawN)/powerLawN))));
				vel[i].x = 0.0f;
				vel[i].y = 0.0;
				pos[i].w = physical_density(vel[i].w, 0)*m_deltap*m_deltap;
			}
		}
	}

void rwN04Fr1Periodic::fillDeviceMap()
	{
		fillDeviceMapByAxis(X_AXIS);
	}


