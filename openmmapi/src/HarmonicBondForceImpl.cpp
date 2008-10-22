/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "internal/OpenMMContextImpl.h"
#include "internal/HarmonicBondForceImpl.h"
#include "kernels.h"

using namespace OpenMM;
using std::pair;
using std::vector;
using std::set;

HarmonicBondForceImpl::HarmonicBondForceImpl(HarmonicBondForce& owner) : owner(owner) {
}

HarmonicBondForceImpl::~HarmonicBondForceImpl() {
}

void HarmonicBondForceImpl::initialize(OpenMMContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcHarmonicBondForceKernel::Name(), context);
    dynamic_cast<CalcHarmonicBondForceKernel&>(kernel.getImpl()).initialize(context.getSystem(), owner);
}

void HarmonicBondForceImpl::calcForces(OpenMMContextImpl& context, Stream& forces) {
    dynamic_cast<CalcHarmonicBondForceKernel&>(kernel.getImpl()).executeForces(context);
}

double HarmonicBondForceImpl::calcEnergy(OpenMMContextImpl& context) {
    return dynamic_cast<CalcHarmonicBondForceKernel&>(kernel.getImpl()).executeEnergy(context);
}

std::vector<std::string> HarmonicBondForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcHarmonicBondForceKernel::Name());
    return names;
}

