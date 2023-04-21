/* Methods of the "ContactDistribution" class.
   Author: Callum Marples
*/

#include "contact_distributions.hpp"
#include "ellipsoid_shape.hpp"
#include "ellipsoid_patches.hpp"
#include "vectors.hpp"
#include <iostream>
#include <vector>
#include <iostream>
#include <memory>

using namespace std;

ContactDistribution::ContactDistribution() {
    cout << "ERROR : REQUIRE ELLIPSOIDSHAPE & ELLIPSOIDPATCHES OBJECTS TO INITIALISE CONTACT DISTRIBUTION" << endl;
    exit(-1);
}

ContactDistribution::ContactDistribution(EllipsoidShape& E, EllipsoidPatches& P) {
    mpShape = make_shared<EllipsoidShape> (E);
    mpPatches = make_shared<PolarPatches> (P);
    mpPatches->InitialiseVector(mContacts);
}

int ContactDistribution::GetNoContacts() {
    int sum = 0;
    for (int i=0; i!=mpPatches->GetNoPatches(); ++i) { sum += mContacts[i]; }
    return sum;
}

void ContactDistribution::AddContact(Vec3& v) {
    // Convert to spherical polars
    double ph = atan2(mpShape->aAxis*v[1], mpShape->bAxis*v[0]);
    if (ph < 0.0) { ph += 2.0*PI; }
    double th = acos(v[2] / mpShape->cAxis);
    UpdateBin(th, ph);
}

void ContactDistribution::UpdateBin(double& th, double& ph) {
    // Find theta and phi indices ...
    int thIndex = mpPatches->FindThetaIndex(th);
    int phIndex = mpPatches->FindPhiIndex(ph);
    // ... and use to find bin index
    int binIndex = mpPatches->PatchIndex(thIndex, phIndex);
    mContacts[binIndex] += 1;
}
