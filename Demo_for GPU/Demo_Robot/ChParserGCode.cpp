// G-Code parser source

#include "ChParserGCode.h"

namespace chrono {

ChParserGCode::ChParserGCode(const std::string& filename, double max_feedrate) {
	m_path = chrono_types::make_shared<ChLinePath>();
	m_feedrates = chrono_types::make_shared<ChFunctionSequence>();
	m_current_frame = m_home_frame;

	m_max_feedrate = max_feedrate;

	try {
		ParseFile(filename);
	}
	catch (const std::exception& exc) {
		std::cerr << "Unable to parse file: " << exc.what() << std::endl;
		throw std::runtime_error("Unable to parse file");
	}
}

std::shared_ptr<ChVisualModel> ChParserGCode::GetVisualModel(const ChFramed& frame) const {
	auto vismodel = chrono_types::make_shared<ChVisualModel>();
	auto path_shape = chrono_types::make_shared<ChVisualShapeLine>();
	path_shape->SetLineGeometry(m_path);
	path_shape->SetColor(ChColor(0.0f, 0.0f, 1.f));
	vismodel->AddShape(path_shape, frame);
	return vismodel;
}

int ChParserGCode::ParseFile(const std::string& filename) {
	std::ifstream filein(filename);
	std::string line = "";
	unsigned int ctr_line = 0;
	// Read file line by line
	while (std::getline(filein, line)) {
		std::vector<std::string> splitted_line = SplitString(line, ' ');
		if (splitted_line.empty()) {
			++ctr_line;
			continue;
		}
		std::string instruction_type = splitted_line[0];

		// Parse single line
		for (const auto& ii : splitted_line) {
			bool continue_line_parsing = ParseSubCommand(ii);
			if (!continue_line_parsing)
				break;
		}
		FinalizeCommand();

		++ctr_line;
	}
	return ctr_line;
}

std::vector<std::string> ChParserGCode::SplitString(const std::string& str, char delimiter) {
	std::stringstream ss(str);
	std::vector<std::string> splitted;
	std::string token;
	while (std::getline(ss, token, delimiter)) {
		if (!token.empty())
			splitted.push_back(token);
	}
	return splitted;
}

double ChParserGCode::GetValueFromSubCommand(std::string subcommand) {
	std::string num_str = subcommand.erase(0, 1);
	return std::stod(num_str);
}

bool ChParserGCode::ParseSubCommand(const std::string& subcommand) {
	switch (subcommand[0]) {
	case ';': // comment
		return false;
		break;
	case '%': // start or end of file
		if (!m_file_started)
			m_file_started = true;
		else
			m_file_ended = true;
		return false;
		break;
	case '(': // comment
		return false;
		break;
	case 'G': { // G subcommand	
		auto m_instr_type_find = m_instructionsetG.find(subcommand);

		if (m_instr_type_find != m_instructionsetG.end())
			m_instr_type = m_instr_type_find->second;
		else {
			std::cerr << "Command: " << subcommand << " not supported." << std::endl;
			return false;
		}

		switch (m_instr_type) {
		case InstructionType::G00:
			m_current_feedrate = length_conversion_factor * m_max_feedrate / 60.0;
			break;
		case InstructionType::G01:
			break;
		case InstructionType::G10:
			break;
		case InstructionType::G17:
			m_working_plane = WorkingPlaneType::XY;
			break;
		case InstructionType::G18:
			m_working_plane = WorkingPlaneType::XZ;
			break;
		case InstructionType::G19:
			m_working_plane = WorkingPlaneType::YZ;
			break;
		case InstructionType::G02:
			break;
		case InstructionType::G20:
			length_conversion_factor = 0.0254;
			break;
		case InstructionType::G21:
			length_conversion_factor = 1e-3;
			break;
		case InstructionType::G28:
			break;
		case InstructionType::G03:
			break;
		case InstructionType::G30:
			break;
		case InstructionType::G38_2:
			break;
		case InstructionType::G04:
			break;
		case InstructionType::G40:
			break;
		case InstructionType::G41:
			break;
		case InstructionType::G42:
			break;
		case InstructionType::G43:
			break;
		case InstructionType::G49:
			break;
		case InstructionType::G53:
			break;
		case InstructionType::G54:
			break;
		case InstructionType::G55:
			break;
		case InstructionType::G56:
			break;
		case InstructionType::G57:
			break;
		case InstructionType::G58:
			break;
		case InstructionType::G59:
			break;
		case InstructionType::G59_1:
			break;
		case InstructionType::G59_2:
			break;
		case InstructionType::G59_3:
			break;
		case InstructionType::G61:
			break;
		case InstructionType::G61_1:
			break;
		case InstructionType::G64:
			break;
		case InstructionType::G80:
			break;
		case InstructionType::G81:
			break;
		case InstructionType::G82:
			break;
		case InstructionType::G83:
			break;
		case InstructionType::G84:
			break;
		case InstructionType::G85:
			break;
		case InstructionType::G86:
			break;
		case InstructionType::G87:
			break;
		case InstructionType::G88:
			break;
		case InstructionType::G89:
			break;
		case InstructionType::G90:
			m_abs_coords = true;
			break;
		case InstructionType::G91:
			m_abs_coords = false;
			break;
		case InstructionType::G92:
			break;
		case InstructionType::G92_1:
			break;
		case InstructionType::G92_2:
			break;
		case InstructionType::G92_3:
			break;
		case InstructionType::G93:
			break;
		case InstructionType::G94:
			break;
		case InstructionType::G98:
			break;
		case InstructionType::G99:
			break;
		default:
			std::cerr << "Developer error: expected a G command." << std::endl;
			break;
		}
		break;
	}
	case 'M': {
		auto m_instr_type_find = m_instructionsetM.find(subcommand);

		if (m_instr_type_find != m_instructionsetM.end())
			m_instr_type = m_instr_type_find->second;
		else {
			std::cerr << "Command: " << subcommand << " not supported." << std::endl;
			return false;
		}

		switch (m_instr_type) {
		case InstructionType::M00:
			break;
		case InstructionType::M01:
			break;
		case InstructionType::M02:
			break;
		case InstructionType::M03:
			break;
		case InstructionType::M30:
			break;
		case InstructionType::M04:
			break;
		case InstructionType::M48:
			break;
		case InstructionType::M49:
			break;
		case InstructionType::M05:
			break;
		case InstructionType::M06:
			break;
		case InstructionType::M60:
			break;
		case InstructionType::M07:
			break;
		case InstructionType::M08:
			break;
		case InstructionType::M09:
			break;
		case InstructionType::M107:
			break;
		default:
			std::cerr << "Developer error: expected a M command." << std::endl;
			break;
		}
		break;
	}
	// Position subcommand
	case 'A':
		break;
	case 'B':
		break;
	case 'C':
		break;
	case 'D':
		break;
	case 'E':
		m_extrusion_length = length_conversion_factor * GetValueFromSubCommand(subcommand);
		break;
	case 'F':
		m_current_feedrate = length_conversion_factor * GetValueFromSubCommand(subcommand) / 60.0;
		break;
	case 'H':
		break;
	case 'I':
		m_arc_center_offset.x() = length_conversion_factor * GetValueFromSubCommand(subcommand);
		break;
	case 'J':
		m_arc_center_offset.y() = length_conversion_factor * GetValueFromSubCommand(subcommand);
		break;
	case 'K':
		m_arc_center_offset.z() = length_conversion_factor * GetValueFromSubCommand(subcommand);
		break;
	case 'L':
		break;
	case 'N':
		break;
	case 'P':
		break;
	case 'Q':
		break;
	case 'R':
		break;
	case 'S':
		m_spindle_speed = GetValueFromSubCommand(subcommand) * CH_RPM_TO_RAD_S;
		break;
	case 'T':
		break;
	case 'X':
		m_abs_coords ? m_next_frame.pos.x() = length_conversion_factor * GetValueFromSubCommand(subcommand) : m_next_frame.pos.x() += length_conversion_factor * GetValueFromSubCommand(subcommand);
		break;
	case 'Y':
		m_abs_coords ? m_next_frame.pos.y() = length_conversion_factor * GetValueFromSubCommand(subcommand) : m_next_frame.pos.y() += length_conversion_factor * GetValueFromSubCommand(subcommand);
		break;
	case 'Z':
		m_abs_coords ? m_next_frame.pos.z() = length_conversion_factor * GetValueFromSubCommand(subcommand) : m_next_frame.pos.z() += length_conversion_factor * GetValueFromSubCommand(subcommand);
		break;
	default:
		std::cerr << "Developer error: command not recognized." << std::endl;
		break;
	}
	return true;
}

void ChParserGCode::FinalizeCommand() {
	auto add_generic_curve = [this]() -> void {
		auto len = m_current_curve->Length(10);
		m_path->AddSubLine(m_current_curve, len);
		auto fun_ramp = chrono_types::make_shared<ChFunctionRamp>(0, m_current_feedrate);
		m_feedrates->InsertFunct(fun_ramp, len / m_current_feedrate, 1.0, true);

		m_current_frame = m_next_frame;
	};

	auto add_arc_curve = [this](bool counterclockwise) {
		ChVector3d arc_center = m_current_frame.pos + m_arc_center_offset;
		ChVector3d arc_frame_X = (m_current_frame.pos - arc_center);

		double radius = arc_frame_X.Length();
		arc_frame_X /= radius;
		ChVector3d arc_dir_Y = (counterclockwise ? 1.0 : -1.0) * (m_next_frame.pos - arc_center).GetNormalized();

		ChVector3d arc_frame_Z = Vcross(arc_frame_X, arc_dir_Y).GetNormalized();
		ChVector3d arc_frame_Y = Vcross(arc_frame_Z, arc_frame_X);
		ChMatrix33d rotmat_glob_loc(arc_frame_X, arc_frame_Y, arc_frame_Z); // from local to global
		ChCoordsysd arc_frame(arc_center, rotmat_glob_loc.GetQuaternion());
		ChVector3d arc_dir_Y_loc = rotmat_glob_loc.transpose() * (m_next_frame.pos - arc_center).GetNormalized();
		double arc_end = atan2(arc_dir_Y_loc.y(), arc_dir_Y_loc.x());
		assert(std::abs(arc_dir_Y_loc.z()) < 1e-5 && "The arc frame computation is ill-conditioned.");
		m_current_curve = chrono_types::make_shared<ChLineArc>(arc_frame, radius, 0, arc_end, counterclockwise);

		double arc_wrap = ChClamp(arc_end, 0., CH_2PI);
		double len = counterclockwise ? radius * (CH_2PI - arc_wrap) : radius * arc_wrap;
		// TODO: move arc length calculation in ChLineArc
		m_path->AddSubLine(m_current_curve, len);
		auto fun_ramp = chrono_types::make_shared<ChFunctionRamp>(0, m_current_feedrate);
		m_feedrates->InsertFunct(fun_ramp, len / m_current_feedrate, 1.0, true);

		m_current_frame = m_next_frame;
	};

	switch (m_instr_type) {
	case InstructionType::G00:
		m_current_curve = chrono_types::make_shared<ChLineSegment>(m_current_frame.pos, m_next_frame.pos);
		add_generic_curve();
		break;
	case InstructionType::G01:
		m_current_curve = chrono_types::make_shared<ChLineSegment>(m_current_frame.pos, m_next_frame.pos);
		add_generic_curve();
		break;
	case InstructionType::G02:
		add_arc_curve(false);
		break;
	case InstructionType::G03:
		add_arc_curve(true);
		break;
	default:
		break;
	}

	m_instr_type = InstructionType::UNKNOWN;
}


} // end namespace chrono