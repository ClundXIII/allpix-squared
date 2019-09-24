/**
 * @file
 * @brief Implements the particle generator
 * @remark Based on code from John Idarraga
 * @copyright Copyright (c) 2017-2020 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GeneratorActionG4.hpp"

#include <limits>
#include <memory>
#include <regex>
#include <sstream>
#include <string>

#include <G4Event.hh>
#include <G4Gamma.hh>
#include <G4IonTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4RunManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4UImanager.hh>
#include <core/module/exceptions.h>

#include "core/config/exceptions.h"
#include "core/utils/log.h"
#include "tools/geant4.h"

using namespace allpix;

std::string allpix::input_energy_file = "";

G4int n_particle = 1;

GeneratorActionG4::GeneratorActionG4(const Configuration& config)
    : fParticleGun(std::make_unique<G4ParticleGun>(n_particle)),
      particle_source_(std::make_unique<G4GeneralParticleSource>()) {

    if(input_energy_file != "") {

        std::cout << "loading stuff from \"" << input_energy_file << "\"" << std::endl;

        std::string base_dir = "../";

        std::ifstream energyInput(base_dir + input_energy_file);

        G4ParticleDefinition* particle = G4Gamma::GammaDefinition();

        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
        fParticleGun->SetParticlePosition(config.get<G4ThreeVector>("source_position"));

        std::string line;

        if(!energyInput.is_open()) {
            std::cout << "cannot open infile!" << std::endl;
        }

        while(std::getline(energyInput, line)) {
            std::vector<double> thisLine;

            std::istringstream iss(line);

            double energy = 0;
            double momX = 0, momY = 0, momZ = 0;

            iss >> energy;
            iss >> momX;
            iss >> momY;
            iss >> momZ;

            thisLine.push_back(energy);
            thisLine.push_back(momX);
            thisLine.push_back(momY);
            thisLine.push_back(momZ);

            particleData.push_back(thisLine);
        }

        energyInput.close();
        std::cout << "done" << std::endl;
        return;
    }

    // Define radioactive isotopes:
    static std::map<std::string, std::tuple<int, int, int, double>> isotopes = {
        {"fe55", std::make_tuple(26, 55, 0, 0.)},
        {"am241", std::make_tuple(95, 241, 0, 0.)},
        {"sr90", std::make_tuple(38, 90, 0, 0.)},
        {"co60", std::make_tuple(27, 60, 0, 0.)},
        {"cs137", std::make_tuple(55, 137, 0, 0.)},
    };

    // Set verbosity of source to off
    particle_source_->SetVerbosity(0);

    // Get source specific parameters
    auto source_type = config.get<std::string>("source_type");

    if(source_type == "macro") {
        LOG(INFO) << "Using user macro for particle source.";

        // Get the UI commander
        G4UImanager* UI = G4UImanager::GetUIpointer();

        // Execute the user's macro
        std::ifstream file(config.getPath("file_name", true));
        std::string line;
        while(std::getline(file, line)) {
            // Check for the "/gps/" pattern in the line:
            if(!line.empty()) {
                if(line.rfind("/gps/number", 0) == 0) {
                    throw ModuleError(
                        "The number of particles must be defined in the main configuration file, not in the macro.");
                } else if(line.rfind("/gps/", 0) == 0 || line.at(0) == '#') {
                    LOG(DEBUG) << "Applying Geant4 macro command: \"" << line << "\"";
                    UI->ApplyCommand(line);
                } else {
                    LOG(WARNING) << "Ignoring Geant4 macro command: \"" + line + "\" - not related to particle source.";
                }
            }
        }

    } else {

        // Get the source and set the centre coordinate of the source
        auto single_source = particle_source_->GetCurrentSource();
        single_source->GetPosDist()->SetCentreCoords(config.get<G4ThreeVector>("source_position"));

        // Set position and direction parameters according to shape
        if(source_type == "beam") {

            // Set position parameters
            single_source->GetPosDist()->SetPosDisType("Beam");
            single_source->GetPosDist()->SetBeamSigmaInR(config.get<double>("beam_size", 0));

            // Set angle distribution parameters
            // NOTE beam2d will always fire in the -z direction of the system
            single_source->GetAngDist()->SetAngDistType("beam2d");

            // Align the -z axis of the system with the direction vector
            auto direction = config.get<G4ThreeVector>("beam_direction");
            if(fabs(direction.mag() - 1.0) > std::numeric_limits<double>::epsilon()) {
                LOG(WARNING) << "Momentum direction is not a unit vector: magnitude is ignored";
            }
            auto min_element = std::min(std::abs(direction.x()), std::min(std::abs(direction.y()), std::abs(direction.z())));
            G4ThreeVector angref1;
            if(min_element == std::abs(direction.x())) {
                angref1 = direction.cross({1, 0, 0});
            } else if(min_element == std::abs(direction.y())) {
                angref1 = direction.cross({0, 1, 0});
            } else if(min_element == std::abs((direction.z()))) {
                angref1 = direction.cross({0, 0, 1});
            }
            G4ThreeVector angref2 = angref1.cross(direction);

            single_source->GetAngDist()->DefineAngRefAxes("angref1", angref1);
            single_source->GetAngDist()->DefineAngRefAxes("angref2", angref2);
            auto divergence = config.get<G4TwoVector>("beam_divergence", G4TwoVector(0., 0.));
            single_source->GetAngDist()->SetBeamSigmaInAngX(divergence.x());
            single_source->GetAngDist()->SetBeamSigmaInAngY(divergence.y());

        } else if(source_type == "sphere") {

            // Set position parameters
            single_source->GetPosDist()->SetPosDisType("Surface");
            single_source->GetPosDist()->SetPosDisShape("Sphere");

            // Set angle distribution parameters
            single_source->GetPosDist()->SetRadius(config.get<double>("sphere_radius"));

            if(config.has("sphere_focus_point")) {
                single_source->GetAngDist()->SetAngDistType("focused");
                single_source->GetAngDist()->SetFocusPoint(config.get<G4ThreeVector>("sphere_focus_point"));
            } else {
                single_source->GetAngDist()->SetAngDistType("cos");
            }

        } else if(source_type == "square") {

            // Set position parameters
            single_source->GetPosDist()->SetPosDisType("Plane");
            single_source->GetPosDist()->SetPosDisShape("Square");
            single_source->GetPosDist()->SetHalfX(config.get<double>("square_side") / 2);
            single_source->GetPosDist()->SetHalfY(config.get<double>("square_side") / 2);

            // Set angle distribution parameters
            single_source->GetAngDist()->SetAngDistType("iso");
            single_source->GetAngDist()->SetMaxTheta(config.get<double>("square_angle", ROOT::Math::Pi()) / 2);

        } else if(source_type == "point") {

            // Set position parameters
            single_source->GetPosDist()->SetPosDisType("Point");

            // Set angle distribution parameters
            single_source->GetAngDist()->SetAngDistType("iso");

        } else {

            throw InvalidValueError(config, "source_type", "");
        }

        // Find Geant4 particle
        auto pdg_table = G4ParticleTable::GetParticleTable();
        auto particle_type = config.get<std::string>("particle_type", "");
        std::transform(particle_type.begin(), particle_type.end(), particle_type.begin(), ::tolower);
        auto particle_code = config.get<int>("particle_code", 0);
        G4ParticleDefinition* particle = nullptr;

        if(!particle_type.empty() && particle_code != 0) {
            if(pdg_table->FindParticle(particle_type) == pdg_table->FindParticle(particle_code)) {
                LOG(WARNING) << "particle_type and particle_code given. Continuing because they match.";
                particle = pdg_table->FindParticle(particle_code);
                if(particle == nullptr) {
                    throw InvalidValueError(config, "particle_code", "particle code does not exist.");
                }
            } else {
                throw InvalidValueError(
                    config, "particle_type", "Given particle_type does not match particle_code. Please remove one of them.");
            }
        } else if(particle_type.empty() && particle_code == 0) {
            throw InvalidValueError(config, "particle_code", "Please set particle_code or particle_type.");
        } else if(particle_code != 0) {
            particle = pdg_table->FindParticle(particle_code);
            if(particle == nullptr) {
                throw InvalidValueError(config, "particle_code", "particle code does not exist.");
            }
        } else if(isotopes.find(particle_type) != isotopes.end()) {
            auto isotope = isotopes[particle_type];
            // Set radioactive isotope:
            particle = G4IonTable::GetIonTable()->GetIon(std::get<0>(isotope), std::get<1>(isotope), std::get<3>(isotope));

            // Force the radioactive isotope to decay immediately:
            particle->SetPDGLifeTime(0.);

            single_source->SetParticleCharge(std::get<2>(isotope));

            // Warn about non-zero source energy:
            if(config.get<double>("source_energy") > 0) {
                LOG(WARNING)
                    << "A radioactive isotope is used as particle source, but the source energy is not set to zero.";
            }
        } else if(particle_type.substr(0, 3) == "ion") {
            // Parse particle type as ion with components /Z/A/Q/E
            std::smatch ion;
            if(std::regex_match(
                   particle_type, ion, std::regex("ion/([0-9]+)/([0-9]+)/([-+]?[0-9]+)/([0-9.]+(?:[a-zA-Z]+)?)")) &&
               ion.ready()) {
                particle = G4IonTable::GetIonTable()->GetIon(
                    allpix::from_string<int>(ion[1]), allpix::from_string<int>(ion[2]), allpix::from_string<double>(ion[4]));
                single_source->SetParticleCharge(allpix::from_string<int>(ion[3]));
            } else {
                throw InvalidValueError(config, "particle_type", "cannot parse parameters for ion.");
            }
        } else {
            particle = pdg_table->FindParticle(particle_type);
            if(particle == nullptr) {
                throw InvalidValueError(config, "particle_type", "particle type does not exist.");
            }
        }

        LOG(DEBUG) << "Using particle " << particle->GetParticleName() << " (ID " << particle->GetPDGEncoding() << ").";

        // Set global parameters of the source
        single_source->SetNumberOfParticles(1);
        single_source->SetParticleDefinition(particle);
        // Set the primary track's start time in for the current event to zero:
        single_source->SetParticleTime(0.0);

        // Set energy parameters
        single_source->GetEneDist()->SetEnergyDisType("Gauss");
        single_source->GetEneDist()->SetMonoEnergy(config.get<double>("source_energy"));
        single_source->GetEneDist()->SetBeamSigmaInE(config.get<double>("source_energy_spread", 0.));
    }
}

/**
 * Called automatically for every event
 */
void GeneratorActionG4::GeneratePrimaries(G4Event* event) {

    if(particleData.size() > 0) {

        fParticleGun->SetParticleEnergy(particleData[particleDataPos][0] * keV);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(
            particleData[particleDataPos][1], particleData[particleDataPos][2], particleData[particleDataPos][3]));
        fParticleGun->GeneratePrimaryVertex(event);

        particleDataPos++;
        if(particleData.size() <= particleDataPos) {
            particleDataPos = 0;
            std::cout << "Input Data wrap around!" << std::endl;
        }
    } else {
        particle_source_->GeneratePrimaryVertex(event);
    }
}
