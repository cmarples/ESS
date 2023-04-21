/* Contact Distribution Class, "ContactDistribution".
   Author: Callum Marples
*/

#ifndef CONTACT_DEF
#define CONTACT_DEF

#include "ellipsoid_shape.hpp"
#include "ellipsoid_patches.hpp"
#include "vectors.hpp"
#include <vector>
#include <memory>

// This class is used to count the number of 'contacts' recorded on each patch
// of an ellipsoid surface.
class ContactDistribution {
    public:
        // Constructors
        ContactDistribution(); // Return error if default used
        ContactDistribution(EllipsoidShape& E, EllipsoidPatches& P);

        // Getters
        int GetNoContacts();
        int operator[](int i) { return mContacts[i]; }
        std::vector<int> GetContactDistribution() { return mContacts; }

        // Add contact
        void AddContact(Vec3& v);
        void AddContact(double& th, double& ph) { UpdateBin(th, ph); }
        void UpdateBin(double& th, double& ph);

    private:
        std::shared_ptr<EllipsoidShape> mpShape;
        std::shared_ptr<PolarPatches> mpPatches;
        std::vector<int> mContacts;
};

#endif // CONTACT_DEF


