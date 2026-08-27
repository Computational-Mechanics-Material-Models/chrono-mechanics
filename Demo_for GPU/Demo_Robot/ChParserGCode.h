// G-Code parser header

#ifndef CH_PARSER_GCODE
#define CH_PARSER_GCODE

#include "chrono/core/ChFrame.h"
#include "chrono/assets/ChVisualModel.h"
#include "chrono/assets/ChVisualShapeLine.h"
#include "chrono/geometry/ChLinePath.h"
#include "chrono/geometry/ChLineArc.h"
#include "chrono/functions/ChFunction.h"
#include "chrono/functions/ChFunctionPositionLine.h"
#include "chrono/utils/ChUtils.h"


namespace chrono {
class ChParserGCode {
public:
	enum class InstructionType {
		A, ///< A-axis of machine
		B, ///< B-axis of machine
		C, ///< C-axis of machine
		D, ///< tool radius compensation number
		E,
		F, ///< feedrate
		H, ///< tool length offset index
		I, ///< X-axis offset for arcs / X offset in G87 canned cycle
		J, ///< Y-axis offset for arcs / Y offset in G87 canned cycle
		K, ///< Z-axis offset for arcs / Z offset in G87 canned cycle
		L, ///< number of repetitions in canned cycles / key used with G10
		N, ///< line number
		P, ///< dwell time in canned cycles / dwell time with G04 / key used with G10
		Q, ///< feed increment in G83 canned cycle
		R, ///< arc radius
		S, ///< canned cycle plane / spindle speed
		T, ///< tool selection
		X, ///< X-axis of machine
		Y, ///< Y-axis of machine
		Z, ///< Z-axis of machine
		G00, ///< rapid positioning
		G01, ///< linear interpolation
		G10, ///< coordinate system origin setting
		G17, ///< XY-plane selection
		G18, ///< XZ-plane selection
		G19, ///< YZ-plane selection
		G02, ///< circular/helical interpolation (clockwise)
		G20, ///< inch system selection
		G21, ///< millimeter system selection
		G28, ///< return to home
		G03, ///< circular/helical interpolation (counterclockwise)
		G30, ///< return to secondary home
		G38_2, ///< straight probe
		G04, ///< dwell
		G40, ///< cancel cutter radius compensation
		G41, ///< start cutter radius compensation left
		G42, ///< start cutter radius compensation right
		G43, ///< tool length offset (plus)
		G49, ///< cancel tool length offset
		G53, ///< motion in machine coordinate system
		G54, ///< use preset work coordinate system 1
		G55, ///< use preset work coordinate system 2
		G56, ///< use preset work coordinate system 3
		G57, ///< use preset work coordinate system 4
		G58, ///< use preset work coordinate system 5
		G59, ///< use preset work coordinate system 6
		G59_1, ///< use preset work coordinate system 7
		G59_2, ///< use preset work coordinate system 8
		G59_3, ///< use preset work coordinate system 9
		G61, ///< set path control mode: exact path
		G61_1, ///< set path control mode: exact stop
		G64, ///< set path control mode: continuous
		G80, ///< cancel motion mode (including any canned cycle)
		G81, ///< canned cycle: drilling
		G82, ///< canned cycle: drilling with dwell
		G83, ///< canned cycle: peck drilling
		G84, ///< canned cycle: right hand tapping
		G85, ///< canned cycle: boring, no dwell, feed out
		G86, ///< canned cycle: boring, spindle stop, rapid out
		G87, ///< canned cycle: back boring
		G88, ///< canned cycle: boring, spindle stop, manual out
		G89, ///< canned cycle: boring, dwell, feed out
		G90, ///< absolute distance mode
		G91, ///< incremental distance mode
		G92, ///< offset coordinate systems and set parameters
		G92_1, ///< cancel offset coordinate systems and set parameters to zero
		G92_2, ///< cancel offset coordinate systems but do not reset parameters
		G92_3, ///< apply parameters to offset coordinate systems
		G93, ///< inverse time feed rate mode
		G94, ///< units per minute feed rate mode
		G98, ///< initial level return in canned cycles
		G99, ///< R-point level return in canned cycles
		M00, ///< program stop
		M01, ///< optional program stop
		M02, ///< program end
		M03, ///< turn spindle clockwise
		M30, ///< program end, pallet shuttle, and reset
		M04, ///< turn spindle counterclockwise
		M48, ///< enable speed and feed overrides
		M49, ///< disable speed and feed overrides
		M05, ///< stop spindle turning
		M06, ///< tool change
		M60, ///< pallet shuttle and program stop
		M07, ///< mist coolant on
		M08, ///< flood coolant on
		M09, ///< mist and flood coolant off
		M107,
		UNKNOWN
	};

	enum class WorkingPlaneType {
		XY, XZ, YZ
	};

	ChParserGCode() = delete;

	/// max_feedrate is assumed to be in the same units as defined in the parsed G-Code
	ChParserGCode(const std::string& filename, double max_feedrate = 100.0);

	std::shared_ptr<ChVisualModel> GetVisualModel(const ChFramed& frame) const;

	int ParseFile(const std::string& filename);

	std::shared_ptr<ChLinePath> GetPath() const { return m_path; }

	std::shared_ptr<ChFunctionSequence> GetSpaceFunction() const { return m_feedrates; }


protected:
	const std::unordered_map<std::string, InstructionType> m_instructionsetG = {
		{"G00", InstructionType::G00},
		{"G0", InstructionType::G00},
		{"G01", InstructionType::G01},
		{"G1", InstructionType::G01},
		{"G10", InstructionType::G10},
		{"G17", InstructionType::G17},
		{"G18", InstructionType::G18},
		{"G19", InstructionType::G19},
		{"G2", InstructionType::G02},
		{"G02", InstructionType::G02},
		{"G20", InstructionType::G20},
		{"G21", InstructionType::G21},
		{"G28", InstructionType::G28},
		{"G03", InstructionType::G03},
		{"G30", InstructionType::G30},
		{"G38.2", InstructionType::G38_2},
		{"G04", InstructionType::G04},
		{"G4", InstructionType::G04},
		{"G40", InstructionType::G40},
		{"G41", InstructionType::G41},
		{"G42", InstructionType::G42},
		{"G43", InstructionType::G43},
		{"G49", InstructionType::G49},
		{"G53", InstructionType::G53},
		{"G54", InstructionType::G54},
		{"G55", InstructionType::G55},
		{"G56", InstructionType::G56},
		{"G57", InstructionType::G57},
		{"G58", InstructionType::G58},
		{"G59", InstructionType::G59},
		{"G59.1", InstructionType::G59_1},
		{"G59.2", InstructionType::G59_2},
		{"G59.3", InstructionType::G59_3},
		{"G61", InstructionType::G61},
		{"G61.1", InstructionType::G61_1},
		{"G64", InstructionType::G64},
		{"G80", InstructionType::G80},
		{"G81", InstructionType::G81},
		{"G82", InstructionType::G82},
		{"G83", InstructionType::G83},
		{"G84", InstructionType::G84},
		{"G85", InstructionType::G85},
		{"G86", InstructionType::G86},
		{"G87", InstructionType::G87},
		{"G88", InstructionType::G88},
		{"G89", InstructionType::G89},
		{"G90", InstructionType::G90},
		{"G91", InstructionType::G91},
		{"G92", InstructionType::G92},
		{"G92.1", InstructionType::G92_1},
		{"G92.2", InstructionType::G92_2},
		{"G92.3", InstructionType::G92_3},
		{"G93", InstructionType::G93},
		{"G94", InstructionType::G94},
		{"G98", InstructionType::G98},
		{"G99", InstructionType::G99},
	};

	const std::unordered_map<std::string, InstructionType> m_instructionsetM = {
		{"M00", InstructionType::M00},
		{"M01", InstructionType::M01},
		{"M02", InstructionType::M02},
		{"M03", InstructionType::M03},
		{"M30", InstructionType::M30},
		{"M04", InstructionType::M04},
		{"M48", InstructionType::M48},
		{"M49", InstructionType::M49},
		{"M05", InstructionType::M05},
		{"M06", InstructionType::M06},
		{"M60", InstructionType::M60},
		{"M07", InstructionType::M07},
		{"M08", InstructionType::M08},
		{"M09", InstructionType::M09},
		{"M107", InstructionType::M107},
		{"M0", InstructionType::M00},
		{"M1", InstructionType::M01},
		{"M2", InstructionType::M02},
		{"M3", InstructionType::M03},
		{"M4", InstructionType::M04},
		{"M5", InstructionType::M05},
		{"M6", InstructionType::M06},
		{"M7", InstructionType::M07},
		{"M8", InstructionType::M08},
		{"M9", InstructionType::M09},
	};

	// Split string separated by given delimiter in vector of substrings
	std::vector<std::string> SplitString(const std::string& str, char delimiter = ' ');

	// Erase first character from a sub-command and return numeric value
	// (e.g. from "X12.3" return 12.3)
	double GetValueFromSubCommand(std::string subcommand);

	bool ParseSubCommand(const std::string& subcommand);

	void FinalizeCommand();

	double length_conversion_factor = 1e-3; ///< conversion factor to apply to parsed numbers to convert to meters
	std::shared_ptr<ChLinePath> m_path = nullptr;
	ChCoordsys<> m_current_frame = CSYSNORM;
	ChCoordsys<> m_next_frame = CSYSNORM;
	ChCoordsys<> m_home_frame = CSYSNORM;
	ChVector3d m_arc_center_offset = VNULL; ///< offset of arc center from starting position
	double m_max_feedrate = 1.0; ///< maximum feedrate (speed) for G00 command; WARNIG: this variable is in [currentunit/min]; must be converted when passed to Chrono
	double m_current_feedrate = 0.0; ///< tool speed [m/s]
	std::shared_ptr<ChLine> m_current_curve = nullptr;
	InstructionType m_instr_type = InstructionType::UNKNOWN;
	bool m_abs_coords = true;
	WorkingPlaneType m_working_plane = WorkingPlaneType::XY;
	bool m_file_started = false;
	bool m_file_ended = false;
	std::shared_ptr<ChFunctionSequence> m_feedrates;
	double m_spindle_speed = 0.0; ///< speed of the spindle [rad/s]
	double m_extrusion_length = 0.0; ///< length of the extruded segment for the current motion [m]
};


} // end namespace chrono

#endif // !CH_PARSER_GCODE
