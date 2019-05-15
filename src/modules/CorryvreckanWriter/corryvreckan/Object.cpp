#include "Object.hpp"
#include "MCParticle.hpp"
#include "Pixel.hpp"

#include "core/utils/exceptions.h"

using namespace corryvreckan;

Object::Object() = default;
Object::Object(std::string detectorID) : m_detectorID(std::move(detectorID)) {}
Object::Object(double timestamp) : m_timestamp(timestamp) {}
Object::Object(std::string detectorID, double timestamp) : m_detectorID(std::move(detectorID)), m_timestamp(timestamp) {}
Object::Object(const Object&) = default;
Object::~Object() = default;

// Return class type for objects which change with detector type
Object* Object::Factory(std::string, std::string objectType, Object* object) {

    if(objectType == "pixels") {
        return (object == nullptr) ? new Pixel() : new Pixel(*static_cast<Pixel*>(object));
    } else if(objectType == "mcparticles") {
        return (object == nullptr) ? new MCParticle() : new MCParticle(*static_cast<MCParticle*>(object));
    }

    return new Object();
}

std::ostream& corryvreckan::operator<<(std::ostream& out, const Object& obj) {
    obj.print(out);
    return out;
}
