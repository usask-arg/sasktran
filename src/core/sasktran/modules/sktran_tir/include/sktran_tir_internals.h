#pragma once
#include <memory>
#include <omp.h>

// SASKTRAN dependencies
#include <modules/sktran_common/sktran_common.h>

// SASKTRAN TIR internal dependencies
#include "sktran_tir_atmosphericstate.h"
#include "sktran_indexofrefraction_tir.h"
#include "sktran_tir_definitions.h"
#include "sktran_tir_raystorage.h"
#include "sktran_tir_rayregionmanager.h"
#include "sktran_tir_opticalproperties.h"
#include "sktran_tir_integrator.h"
#include "sktran_tir_wf.h"
#include "sktran_tir_specs_user_integrator.h"
#include "sktran_tir_specs_user_raytracer.h"
#include "sktran_tir_specs_user_opticalpropertiestable.h"
#include "sktran_tir_specs_user_wf.h"
#include "sktran_tir_specs_user.h"
#include "sktran_tir_linesofsighttable.h"
#include "sktran_tir_specs_internal_raytracer.h"
#include "sktran_tir_specs_internal_opticalpropertiestable.h"
#include "sktran_tir_specs_internal_integrator.h"
#include "sktran_tir_specs_internal_wf.h"
#include "sktran_tir_specs_internal_core.h"
#include "sktran_tir_thread_manager.h"
#include "sktran_tir_engine.h"
