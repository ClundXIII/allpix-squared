/**
 * @file
 * @brief Implementation of module to read doping concentration maps
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "DopingProfileReaderModule.hpp"

#include <string>
#include <utility>

#include "core/utils/log.h"

using namespace allpix;

DopingProfileReaderModule::DopingProfileReaderModule(Configuration& config, Messenger*, std::shared_ptr<Detector> detector)
    : Module(config, detector), detector_(std::move(detector)) {}

void DopingProfileReaderModule::init() {
    FieldType type = FieldType::GRID;

    // Check field strength
    auto field_model = config_.get<std::string>("model");

    // Calculate the field depending on the configuration
    if(field_model == "init" || field_model == "apf") {
        // Read the field scales from the configuration, defaulting to 1.0x1.0 pixel cell:
        auto scales = config_.get<ROOT::Math::XYVector>("field_scale", {1.0, 1.0});
        // FIXME Add sanity checks for scales here
        LOG(DEBUG) << "Doping concentration map will be scaled with factors " << scales;
        std::array<double, 2> field_scale{{scales.x(), scales.y()}};

        // Get the field offset in fractions of the pixel pitch, default is 0.0x0.0, i.e. starting at pixel boundary:
        auto offset = config_.get<ROOT::Math::XYVector>("field_offset", {0.0, 0.0});
        if(offset.x() > 1.0 || offset.y() > 1.0) {
            throw InvalidValueError(
                config_,
                "field_offset",
                "shifting doping concentration map by more than one pixel (offset > 1.0) is not allowed");
        }
        LOG(DEBUG) << "Doping concentration map starts with offset " << offset << " to pixel boundary";
        std::array<double, 2> field_offset{{offset.x(), offset.y()}};

        auto field_data = read_field(field_scale, field_model);
        detector_->setDopingProfileGrid(field_data.getData(), field_data.getDimensions(), field_scale, field_offset);

    } else if(field_model == "constant") {
        LOG(TRACE) << "Adding constant doping concentration";
        type = FieldType::CONSTANT;

        auto concentration = config_.get<double>("doping_concentration");
        LOG(INFO) << "Set constant doping concentration of " << Units::display(concentration, {"/cm/cm/cm"});
        FieldFunction<double> function = [concentration](const ROOT::Math::XYZPoint&) { return concentration; };

        detector_->setDopingProfileFunction(function, type);
    } else {
        throw InvalidValueError(config_, "model", "model should be 'constant' or 'init'");
    }
}

/**
 * The field read from the INIT format are shared between module instantiations using the static FieldParser.
 */
FieldParser<double> DopingProfileReaderModule::field_parser_(FieldQuantity::SCALAR);
FieldData<double> DopingProfileReaderModule::read_field(std::array<double, 2> field_scale, const std::string& format) {

    FileType type = (format == "init" ? FileType::INIT : format == "apf" ? FileType::APF : FileType::UNKNOWN);
    std::string units = (type == FileType::INIT ? "/cm/cm/cm" : "");

    try {
        LOG(TRACE) << "Fetching doping concentration map from init file";

        // Get field from file
        auto field_data = field_parser_.get_by_file_name(config_.getPath("file_name", true), type, units);

        // Check if electric field matches chip
        check_detector_match(field_data.getSize(), field_scale);

        LOG(INFO) << "Set doping concentration map with " << field_data.getDimensions().at(0) << "x"
                  << field_data.getDimensions().at(1) << "x" << field_data.getDimensions().at(2) << " cells";

        // Return the field data
        return field_data;
    } catch(std::invalid_argument& e) {
        throw InvalidValueError(config_, "file_name", e.what());
    } catch(std::runtime_error& e) {
        throw InvalidValueError(config_, "file_name", e.what());
    } catch(std::bad_alloc& e) {
        throw InvalidValueError(config_, "file_name", "file too large");
    }
}

/**
 * @brief Check if the detector matches the file header
 */
void DopingProfileReaderModule::check_detector_match(std::array<double, 3> dimensions, std::array<double, 2> field_scale) {
    auto xpixsz = dimensions[0];
    auto ypixsz = dimensions[1];
    auto thickness = dimensions[2];

    auto model = detector_->getModel();
    // Do a several checks with the detector model
    if(model != nullptr) {
        // Check field dimension in z versus the sensor thickness:
        if(std::fabs(thickness - model->getSensorSize().z()) > std::numeric_limits<double>::epsilon()) {
            LOG(WARNING) << "Thickness of doping concentration map is " << Units::display(thickness, "um")
                         << " but sensor thickness is " << Units::display(model->getSensorSize().z(), "um");
        }

        // Check the field extent along the pixel pitch in x and y:
        if(std::fabs(xpixsz - model->getPixelSize().x() * field_scale[0]) > std::numeric_limits<double>::epsilon() ||
           std::fabs(ypixsz - model->getPixelSize().y() * field_scale[1]) > std::numeric_limits<double>::epsilon()) {
            LOG(WARNING) << "Doping concentration map size is (" << Units::display(xpixsz, {"um", "mm"}) << ","
                         << Units::display(ypixsz, {"um", "mm"}) << ") but current configuration results in an map area of ("
                         << Units::display(model->getPixelSize().x() * field_scale[0], {"um", "mm"}) << ","
                         << Units::display(model->getPixelSize().y() * field_scale[1], {"um", "mm"}) << ")" << std::endl
                         << "The size of the area to which the doping concentration is applied can be changes using the "
                            "field_scale parameter.";
        }
    }
}