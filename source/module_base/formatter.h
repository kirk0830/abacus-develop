#ifndef FORMATTER_H
#define FORMATTER_H

#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <utility>
/*
* # ABACUS formatter library
* providing from basic the single variable formatting to complex vectors, tables formatting output
* Author: Yike HUANG, Daye ZHENG
* Institution: AISI
* ps. Zuxin JIN also helps solve dependency problems
* ## Classes
* - Fmt: single variable formatting (0 dimension)
* - PhysicalFmt: bundle with Fmt object, work as a decorator, adjusting the Fmt object according to the physical context, energy, force, coordinate, etc.
* - Table: output data structure, can be used to establish a table with titles, frames, etc.
* - ContextFmt: context formatting, providing the most sophisticated formatting
* ## Usage
* There is a Chinese Feishu document available for this library, please refer to it for more details:
* https://ucoyxk075n.feishu.cn/docx/Yym9dnm3aoTMfHxin8rcX9Rvnmb?theme=LIGHT&contentTheme=DARK
* ## Structure
* - File structure:
* formatter.h: header file
* formatter.cpp: source file
* test/formatter_fmt_test.cpp: unit test for Fmt class
* test/formatter_physfmt_test.cpp: unit test for PhysicalFmt class
* test/formatter_table_test.cpp: unit test for Table class
* test/formatter_contextfmt_test.cpp: unit test for ContextFmt class
* - Class structure:
*               Fmt: single variable formatting (0 dimension)
*                |
*                |
*                |
*          PhysicalFmt: decorator of Fmt
*                |
*                |
*                |
*           ContextFmt <----(inheritance)---- Table
*/
namespace formatter
{
    class Fmt {
        public:
            /// @brief default constructor, set default values: width = 4, precision = 2, fillChar = ' ', fixed = true, right = true, error = false
            Fmt();
            /// @brief constructor with input values
            /// @param width width of the output string
            /// @param precision precision of the output string
            /// @param fillChar fill character of the output string
            /// @param fixed whether the output string keeps decimal places
            /// @param right whether the output string is right aligned
            /// @param error whether the output string has "+" sign if positive
            Fmt(int width, int precision, char fillChar, bool fixed, bool right, bool error);
            /// @brief destructor
            ~Fmt();
            /// @brief setter the width
            /// @param width 
            void set_width(int width) { width_ = width; }
            /// @brief setter the precision
            /// @param precision
            void set_precision(int precision) { precision_ = precision; }
            /// @brief setter the fillChar
            /// @param fillChar
            void set_fillChar(char fillChar) { fillChar_ = fillChar; }
            /// @brief setter the fixed
            /// @param fixed
            void set_fixed(bool fixed) { fixed_ = fixed; }
            /// @brief setter the right
            /// @param right
            void set_right(bool right) { right_ = right; }
            /// @brief setter the error
            /// @param error
            void set_error(bool error) { error_ = error; }

            /// @brief format the input value to string
            /// @tparam T double, int or std::string
            /// @param value input value
            /// @return std::string, the formatted string
            template <typename T> std::string format(const T& value);
            /// @brief getter the width
            /// @return width
            int get_width() const { return width_; }
            /// @brief getter the precision
            /// @return precision
            int get_precision() const { return precision_; }
            /// @brief getter the fillChar
            /// @return fillChar
            char get_fillChar() const { return fillChar_; }
            /// @brief getter the fixed
            /// @return fixed
            bool get_fixed() const { return fixed_; }
            /// @brief getter the right
            /// @return right
            bool get_right() const { return right_; }
            /// @brief getter the error
            /// @return error
            bool get_error() const { return error_; }

        private:
            /// @brief width of the output string
            int width_ = 4;
            /// @brief precision of the output string
            int precision_ = 2;
            /// @brief fill character of the output string
            char fillChar_ = ' ';
            /// @brief whether the output string keeps decimal places
            bool fixed_ = true;
            /// @brief whether the output string is right aligned
            bool right_ = true;
            /// @brief whether the output string has "+" sign if positive
            bool error_ = false;
    };

    class PhysicalFmt {
        public:
            /// @brief default constructor, set default values: context = "none", p_formatter = nullptr
            /// @param context available choices see function adjust_formatter()
            /// @param p_formatter pointer to the formatter, can be passed with pointer to Fmt object or nullptr
            PhysicalFmt(std::string context = "none", Fmt* p_formatter = nullptr);
            /// @brief destructor
            ~PhysicalFmt();

            /// @brief adjust the formatter according to the context
            /// @param left whether the output string is left aligned
            void adjust_formatter(bool left = false);
            /// @brief setter the context
            /// @param context available choices see function adjust_formatter()
            void set_context(std::string context);
            /// @brief setter the pointer to the formatter
            /// @param p_formatter pointer to the formatter, can be passed with pointer to Fmt object or nullptr
            void set_p_formatter(Fmt* p_formatter);
            /// @brief setter the decorator mode, only for debug and unittest purpose
            /// @param decorator_mode if external Fmt object is set to bundle with PhysicalFmt object, then decorator_mode should be set true, otherwise false
            void set_decorator_mode(bool decorator_mode) { decorator_mode_ = decorator_mode; }
            /// @brief getter the context
            /// @return context
            std::string get_context() const { return context_; }
            /// @brief getter the pointer to the formatter
            /// @return pointer to the formatter
            Fmt* get_p_formatter() const { return p_formatter_; }
            /// @brief getter the decorator mode, only for debug and unittest purpose
            /// @return decorator_mode
            bool get_decorator_mode() const { return decorator_mode_; }
        private:
            /// @brief context, available choices see function adjust_formatter()
            std::string context_ = "none";
            /// @brief pointer to the formatter, can be passed with pointer to Fmt object or nullptr
            Fmt* p_formatter_;

            /// @brief decorator mode, indicating whether the external Fmt object is set to bundle with PhysicalFmt object
            bool decorator_mode_ = false;
    };

    class Table 
    {
        public:
            /// @brief default constructor
            Table(){};
            /// @brief destructor
            ~Table(){};
            /// @brief setter the mode
            /// @param new_mode 0: table, 1: data, -1: header
            void set_mode(int new_mode) { mode_ = new_mode; }
            /// @brief setter the overall title
            /// @param new_title overall title
            void set_overall_title(std::string new_title) { overall_title = new_title; }
            /// @brief setter the titles of each column
            /// @param new_titles std::vector<std::string> of titles
            void set_titles(std::vector<std::string> new_titles) { titles_ = new_titles; }
            /// @brief setter the title for one specific column
            /// @param icol column index
            /// @param new_title new title
            void set_title(int icol, std::string new_title) { titles_[icol] = new_title; }
            /// @brief setter the column delimiter
            /// @param new_delimiter new delimiter
            void set_col_delimiter(char new_delimiter) { col_delimiter_ = new_delimiter; }
            /// @brief setter the frame switches
            /// @param new_switches std::vector<int> of switches, 1: print frame, 0: not print frame
            void set_frame_switches(std::vector<int> new_switches) { frame_switches_ = new_switches; }
            /// @brief setter the frame switch for one specific frame
            /// @param frame frame index, 0: top, 1: down, 2: left, 3: right
            /// @param new_switch new switch, 1: print frame, 0: not print frame
            void set_frame_switch(int frame, int new_switch) { frame_switches_[frame] = new_switch; }
            /// @brief setter the frame delimiters
            /// @param new_delimiters std::vector<char> of delimiters
            void set_frame_delimiters(std::vector<char> new_delimiters) { frame_delimiters_ = new_delimiters; }
            /// @brief setter the frame delimiter for one specific frame
            /// @param iframe frame index, 0: top, 1: down, 2: left, 3: right
            /// @param new_delimiter new delimiter
            void set_frame_delimiter(int iframe, char new_delimiter) { frame_delimiters_[iframe] = new_delimiter; }
            /// @brief setter the frame mid switch
            /// @param new_switch new switch, 1: print frame, 0: not print frame
            void set_frame_mid_switch(int new_switch) { frame_mid_switch_ = new_switch; }
            /// @brief setter the frame mid delimiter
            /// @param new_delimiter new delimiter
            void set_frame_mid_delimiter(char new_delimiter) { frame_mid_delimiter_ = new_delimiter; }

            /// @brief enable the flexible width, if set true, for each column different width will be used
            void enable_flexible_width() { flexible_width_ = true; }
            /// @brief disable the flexible width
            void disable_flexible_width() { flexible_width_ = false; }
            /// @brief print frame for top
            void enable_up_frame() { frame_switches_[0] = 1; }
            /// @brief not print frame for top
            void disable_up_frame() { frame_switches_[0] = 0; }
            /// @brief print frame for down
            void enable_down_frame() { frame_switches_[1] = 1; }
            /// @brief not print frame for down
            void disable_down_frame() { frame_switches_[1] = 0; }
            /// @brief print frame for left
            void enable_left_frame() { frame_switches_[2] = 1; }
            /// @brief not print frame for left
            void disable_left_frame() { frame_switches_[2] = 0; }
            /// @brief print frame for right
            void enable_right_frame() { frame_switches_[3] = 1; }
            /// @brief not print frame for right
            void disable_right_frame() { frame_switches_[3] = 0; }
            /// @brief print frame between title and data
            void enable_mid_frame() { frame_mid_switch_ = 1; }
            /// @brief not print frame between title and data
            void disable_mid_frame() { frame_mid_switch_ = 0; }

            /// @brief centerize the title
            void center_title() { title_centered_ = true; }
            /// @brief add new column of data to the table
            /// @param new_title title of the new column
            /// @param new_col data, stored in std::vector<std::string>
            /// @return return the index of the new column
            int add_col(std::string new_title, std::vector<std::string> new_col);
            int add_cols(std::vector<std::string> new_titles, std::vector<std::vector<std::string>> new_cols); // later
            void del_col(int col);  // not implemented
            void del_col(std::string title); // not implemented
            void del_cols(std::vector<std::string> titles); // not implemented

            /// @brief adjust the width of each column according to the data and title
            void adjust_col_width();
            /// @brief DO NOT CALL IT, it is called by adjust_col_width()
            void centerize_title();
            /// @brief clean the data
            void clean();
            /// @brief reset the table
            void reset();
            /// @brief print the table
            std::string print_table();

            /// @brief getter the mode
            /// @return mode, 0: table, 1: data, -1: header
            int get_mode() const { return mode_; }
            /// @brief getter the overall title
            /// @return overall title
            std::string get_overall_title() const { return overall_title; }
            /// @brief getter the number of columns
            /// @return number of columns
            int get_ncol() const { return ncol_; } 
            /// @brief getter the titles of each column
            /// @return std::vector<std::string> of titles
            std::vector<std::string> get_titles() const { return titles_; }
            /// @brief getter the delimiter between columns
            /// @return delimiter
            char get_col_delimiter() const { return col_delimiter_; }
            /// @brief getter the frame switches
            /// @return std::vector<int> of switches, 1: print frame, 0: not print frame
            std::vector<int> get_frame_switches() const { return frame_switches_; }
            /// @brief getter the frame delimiters
            /// @return std::vector<char> of delimiters
            std::vector<char> get_frame_delimiters() const { return frame_delimiters_; }
            /// @brief getter the frame mid switch
            /// @return frame mid switch, 1: print frame, 0: not print frame
            int get_frame_mid_switch() const { return frame_mid_switch_; }
            /// @brief getter the frame mid delimiter
            /// @return frame mid delimiter
            char get_frame_mid_delimiter() const { return frame_mid_delimiter_; }
            /// @brief getter the boolean controlling flexible width switch
            /// @return flexible width switch
            bool get_flexible_width() const {return flexible_width_;}
            /// @brief getter the width of each column
            /// @return std::vector<int> of width
            std::vector<int> get_col_widths() const { return col_widths_; }
            /// @brief getter the maximum width of each column
            /// @return maximum width
            int get_col_max_width() const { return col_max_width_; }
            /// @brief getter the total width of the table, will be the sum of column width and delimiter width plus 3 (delimiters)
            /// @return total width
            int get_total_width() const { return total_width_; }
            /// @brief getter the data of the table
            /// @return std::vector<std::vector<std::string>> of data
            std::vector<std::vector<std::string>> get_data() const { return data_; }
        private:
            int mode_ = 0; // 0: table, 1: data, -1: header

            /// @brief overall title
            std::string overall_title;

            /// @brief number of columns
            int ncol_ = 0;
            /// @brief titles of each column
            std::vector<std::string> titles_;
            /// @brief delimiter between columns
            char col_delimiter_ = ' ';
            /// @brief frame switches, 1: print frame, 0: not print frame, positional info.: 0: top, 1: down, 2: left, 3: right
            std::vector<int> frame_switches_ = {1, 1, 0, 0};
            /// @brief frame delimiters, positional info.: 0: top, 1: down, 2: left, 3: right
            std::vector<char> frame_delimiters_ = {'-', '-', '|', '|'};
            /// @brief frame mid switch, the mid frame means the frame between title and data, 1: print frame, 0: not print frame
            int frame_mid_switch_ = 1;
            /// @brief frame mid delimiter
            char frame_mid_delimiter_ = '-';
            /// @brief if set true, for each column different width will be used
            bool flexible_width_ = true;
            /// @brief width of each column
            std::vector<int> col_widths_;
            /// @brief maximum width of each column
            int col_max_width_ = 0;
            /// @brief total width of the table
            int total_width_ = 0;
            /// @brief data of the table
            bool title_centered_ = false;
            // warning: the col data is stored here row-by-row
            std::vector<std::vector<std::string>> data_;
    };

    class ContextFmt : public Table
    {
        public:
            /// @brief default constructor, for default values see declarations
            ContextFmt();
            /// @brief destructor
            ~ContextFmt();

            // physical formatter default context: energy
            /// @brief set context by context name, if it is already a predefined type, then directly initialize phys_fmt_ vector, set number of columns and title.
            /// @param context context name, see this->predefined_phys_fmt
            /// @note this function is mostly used for inputing data in 0 dimension one-by-one
            /// @attention for pointer pass to << operator, you MUST NOT USE THIS but to provide another int* nrows!
            void set_context(std::string context); // most simple case, rely on default values pre-defined
            /// @brief (pointer overloaded version) set context by context name, , if it is already a predefined type, then directly initialize phys_fmt_ vector, set number of columns and title.
            /// @param context context name, see this->predefined_phys_fmt
            /// @param ncol number of columns, needed to tranverse the last param nrows
            /// @param nrows (pointer case specifically compulsory) rows of every column
            void set_context(std::string context, int ncol, int* &nrows);
            /// @brief (pointer overloaded version) set context directly by properties of each column stored in std::string*, also set number of rows for each column
            /// @param ncol number of columns
            /// @param phys_fmt physical format of each column
            /// @param nrows number of rows for each column
            void set_context(int ncol, std::string* phys_fmt, int* nrows); // general case, for array input
            /// @brief (STL vector overloaded version) set context directly by properties of each column stored in std::vector<std::string>, the number of rows for each column is automatically contained in input data vector.
            /// @param phys_fmt physical format of each column
            void set_context(std::vector<std::string> phys_fmt); // general case, for vector input
            
            /// @brief setter the default physical context, for cases the context is not pre defined, like vector3d, vector3d_i, etc.
            /// @param default_phys_fmt default physical context, default is "energy" now
            void set_default_phys_fmt(std::string default_phys_fmt) { default_phys_fmt_ = default_phys_fmt; }
            /// @brief setter the cache title, to be used for adding new column to table
            /// @param cache_title cache title
            void set_cache_title(std::string cache_title) { cache_title_ = cache_title; }
            /// @brief setter the number of columns, only used for pointer input data
            /// @param ncol number of columns
            void set_ncol(int ncol) { ncol_ = ncol; }
            /// @brief setter the current column index
            /// @param icol current column index
            void set_icol(int icol) { icol_ = icol; }

            /// @brief setter the title_switch
            /// @param title_switch title switch, 0: not print title, 1: print title
            void set_title_switch(int title_switch) { title_switch_ = title_switch; }

            /// @brief overloaded operator<<, for inputing data in 0 dimension one-by-one
            /// @tparam T int, double, std::string
            /// @param value input value
            /// @return *this
            /// @attention if enable_title(), you can input title and then the data one-by-one, but if disable_title(), you can only input data one-by-one, or unexpected result will be produced
            template <typename T>
            formatter::ContextFmt& operator<<(const T& value);
            /// @brief overloaded operator<<, for inputing data in 1 dimension std::vector<T>
            /// @tparam T int, double and std::string
            /// @param value input value
            /// @return *this
            /// @attention if enable_title(), you can input title and then the data one-by-one, but if disable_title(), you can only input data one-by-one, or unexpected result will be produced
            template <typename T>
            formatter::ContextFmt& operator<<(const std::vector<T>& value);
            /// @brief overloaded operator<<, for inputing data in 1 dimension T*
            /// @tparam T int, double and std::string
            /// @param value input value
            /// @return *this
            /// @attention if enable_title(), you can input title and then the data one-by-one, but if disable_title(), you can only input data one-by-one, or unexpected result will be produced
            template <typename T>
            formatter::ContextFmt& operator<<(T*& value);

            /// @brief overloaded operator<<, for inputing data in char*
            /// @param value input value
            /// @return *this
            formatter::ContextFmt& operator<<(char const* value);

            /// @brief reset all
            void reset();
            /// @brief print the table
            std::string str(bool next_line = false);
            /// @brief set to mode in which title will be manually input and certainly will be output
            void enable_title() { Table::set_mode(0); this->set_title_switch(0); }
            /// @brief set to mode in which title will not be manually input and certainly will not be output
            void disable_title() { Table::set_mode(1); this->set_title_switch(1); }
            /// @brief actually it is something strange and the function is not recommended to use
            void only_title() { Table::set_mode(-1); this->set_title_switch(0); }
            /// @brief enable the iterative mode, in which the data will be output in iterative way, for details see online document whose link is provided in the header of this file
            void enable_iterative() { iterative_ = 1; }
            /// @brief disable the iterative mode
            void disable_iterative() { iterative_ = 0; }

            /// @brief getter the context
            /// @return context
            std::string get_context() const { return context_; }
            /// @brief getter the physical format of each column
            /// @return std::vector<std::string> of physical format
            std::vector<std::string> get_phys_fmt() const { return phys_fmt_; }
            /// @brief getter the pointer to the physical formatter
            /// @return pointer to the physical formatter
            formatter::PhysicalFmt* get_p_phys_fmt() { return p_phys_fmt_; }
            /// @brief getter the default physical context
            /// @return default physical context
            std::string get_default_phys_fmt() const { return default_phys_fmt_; }
            /// @brief getter the cache title
            /// @return cache title
            std::string get_cache_title() const { return cache_title_; }
            /// @brief getter the number of columns, only used for pointer input data
            /// @return number of columns
            int get_ncol() const { return ncol_; }
            /// @brief getter the current column index
            /// @return current column index
            int get_icol() const { return icol_; }
            /// @brief getter the number of rows for each column
            /// @return std::vector<int> of number of rows
            std::vector<int> get_nrows() const { return nrows_; }
            /// @brief getter the formatter
            /// @return formatter
            formatter::Fmt& get_fmt() { return fmt_; }
            /// @brief getter the stringstream (did I use it? I don't remember)
            /// @return stringstream
            std::stringstream& get_ss() { return ss_; }
            /// @brief getter the title switch
            /// @return title switch
            int get_title_switch() const { return title_switch_; }
            /// @brief getter the iterative mode
            /// @return iterative mode
            bool get_iterative() const { return iterative_; }
            // for debug, no need to test
            void print_status() const;
        private:
            /// @brief context
            std::string context_;
            /// @brief physical format of each column
            std::vector<std::string> phys_fmt_;
            /// @brief default physical context
            std::string default_phys_fmt_ = "energy";
            /// @brief cache title, to be used for adding new column to table
            std::string cache_title_ = "";
            /// @brief whether the title is manually input
            bool with_title_ = true;

            /// @brief number of columns, only used for pointer input data
            int ncol_ = 0;
            /// @brief current column index
            int icol_ = 0;
            /// @brief number of rows for each column
            std::vector<int> nrows_; // for operator<< the T* input case
            /// @brief boolean controlling whether the iterative mode is enabled
            int iterative_ = 0;

            /// @brief formatter
            formatter::Fmt fmt_;
            /// @brief pointer to the physical formatter
            formatter::PhysicalFmt* p_phys_fmt_;
            /// @brief stringstream
            std::stringstream ss_;
            /// @brief title switch, 0: not print title, 1: print title
            int title_switch_ = 0;
            // values of title_switch_:
            // title: 0 2 4 6 8 ... , for mode = 1, stepsize = 2, otherwise 1
            // data:  1 3 5 7 9 ...

            /// @brief boolean controlling whether the context is known
            bool known_context_ = false;
            /// @brief database of predefined physical format
            std::unordered_map<std::string, std::pair<int, std::vector<std::string>>> predefined_phys_fmt = {
                {"vector3d", std::make_pair(1,
                std::vector<std::string>{"coordinate", "coordinate", "coordinate"})}, // one can define general context (recommended)
                {"vector3d_i", std::make_pair(1,
                std::vector<std::string>{"constraint", "constraint", "constraint"})}, // vector3d will be position, vectors and for this, it is constraint, kmesh, ...
                {"scf", std::make_pair(0,
                std::vector<std::string>{"str_w4", "int_w4", "energy", "energy", "energy", "time"})}, // but for scf it is really a special case
                {"time_statistics", std::make_pair(0,
                std::vector<std::string>{"str_w30", "str_w30", "time", "int_w8", "double_w6_f2", "double_w6_f2"})}, // so is time statistics
                {"atomic_species", std::make_pair(1,
                std::vector<std::string>{"str_w4", "mass", "str_w30"})}, // so is ATOMIC_SPECIES
                {"lattice_vectors", std::make_pair(1,
                std::vector<std::string>{"coordinate", "coordinate", "coordinate"})} // but it is not true for LATTICE_VECTORS
            };
    };
}
/// @brief global variable, for convenience, the same as INPUT
extern formatter::ContextFmt context;
#endif
