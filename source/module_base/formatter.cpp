#include "formatter.h"
#include <cstring>

formatter::Fmt::Fmt() {
    width_ = 4;
    precision_ = 2;
    fillChar_ = ' ';
    fixed_ = true;
    right_ = true;
    error_ = false;
}

formatter::Fmt::Fmt(int width, int precision, char fillChar, bool fixed, bool right, bool error) {
    width_ = width;
    precision_ = precision;
    fillChar_ = fillChar;
    fixed_ = fixed;
    right_ = right;
    error_ = error;
}

formatter::Fmt::~Fmt() {}

template <typename T>
std::string formatter::Fmt::format(const T& value) {
    std::stringstream ss_;
    if (std::is_signed<T>::value) {
        if (value >= 0) {
            if (error_) ss_ << "+";
        } else {}
    }
    ss_ << std::setprecision(precision_);
    std::stringstream ss;
    if (fixed_) {
        ss_ << std::fixed;
    } else {
        ss_ << std::scientific;
    }
    ss_ << (double)value;

    if (!right_) ss << std::left;
    ss << std::setw(width_) << std::setfill(fillChar_);
    ss << ss_.str();
    return ss.str();
}

template<>
std::string formatter::Fmt::format(const std::string& value) {
    std::stringstream ss;
    if (!right_) ss << std::left;
    ss << std::setw(width_) << std::setfill(fillChar_);
    ss << value;
    return ss.str();
}

// it's not good practice to overload such a function in a relative basic class, which will spoil the whole idea of inheritance
// so here I comment them two out
/*
template <typename T>
std::string formatter::Fmt::format(const std::vector<T>& value) {
    std::stringstream ss;
    for (auto v : value) {
        ss << this->format(v);
    }
    return ss.str();
}

template <typename T>
std::string formatter::Fmt::format(const T* value, int size) {
    std::stringstream ss;
    for (int i = 0; i < size; i++) {
        ss << this->format(value[i]);
    }
    return ss.str();
}
*/

formatter::PhysicalFmt::PhysicalFmt(std::string context, Fmt* p_formatter) {
    context_ = context;
    if (p_formatter != nullptr) {
        this->p_formatter_ = p_formatter;
        this->decorator_mode_ = true;
    }
    else {
        this->p_formatter_ = new Fmt();
    }
    this->adjust_formatter();
}

formatter::PhysicalFmt::~PhysicalFmt() {

    if (this->p_formatter_ != nullptr && !this->decorator_mode_) {
        delete this->p_formatter_;
    }
}

void formatter::PhysicalFmt::adjust_formatter(bool left) {

    auto context = this->context_.c_str();
    if (strcmp(context, "none") == 0) {
        return;
    }
    else if (
        (strcmp(context, "int_w2") == 0)
      ||(strcmp(context, "kmesh") == 0)
      ||(strcmp(context, "constraint") == 0)
    ) {
        this->p_formatter_->set_width(2); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "int_w4") == 0)
      ||(strcmp(context, "scf_step") == 0)
    ) {
        this->p_formatter_->set_width(4); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(left);
    }
    else if (
        (strcmp(context, "int_w8") == 0)
    ) {
        this->p_formatter_->set_width(8); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w6_f1") == 0)
      ||(strcmp(context, "time") == 0)
        ) {
        this->p_formatter_->set_width(6); this->p_formatter_->set_precision(1);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w6_f2") == 0)
      ||(strcmp(context, "mass") == 0)
    ) {
        this->p_formatter_->set_width(6); this->p_formatter_->set_precision(2);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w8_f4_scientific") == 0)
      ||(strcmp(context, "threshold") == 0)
    ) {
        this->p_formatter_->set_width(8); this->p_formatter_->set_precision(4);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(false);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w8_f4_error") == 0)
      ||(strcmp(context, "charge") == 0)
    ) {
        this->p_formatter_->set_width(6); this->p_formatter_->set_precision(4);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left); this->p_formatter_->set_error(true);
    }
    else if (
        (strcmp(context, "double_w10_f2") == 0)
    ) {
        this->p_formatter_->set_width(10); this->p_formatter_->set_precision(2);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w16_f10") == 0)
      ||(strcmp(context, "coordinate") == 0)
      ||(strcmp(context, "position") == 0)
      ||(strcmp(context, "displacement") == 0)
      ||(strcmp(context, "lattice_constant") == 0)
      ||(strcmp(context, "lattice_vector") == 0)
      ||(strcmp(context, "lattice") == 0)
      ||(strcmp(context, "velocity") == 0)
    ) {
        this->p_formatter_->set_width(16); this->p_formatter_->set_precision(10);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "double_w20_f10") == 0)
      ||(strcmp(context, "force") == 0)
      ||(strcmp(context, "energy") == 0)
    ) {
        this->p_formatter_->set_width(20); this->p_formatter_->set_precision(10);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "str_w4") == 0)
      ||(strcmp(context, "numbered_item") == 0)
    ) {
        this->p_formatter_->set_width(4); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(!left);
    }
    else if (
        (strcmp(context, "str_w30") == 0)
      ||(strcmp(context, "long_title") == 0)
    ) {
        this->p_formatter_->set_width(30); this->p_formatter_->set_precision(0);
        this->p_formatter_->set_fillChar(' '); this->p_formatter_->set_fixed(true);
        this->p_formatter_->set_right(left);
    }
    else {
    }
}

void formatter::PhysicalFmt::set_context(std::string context) {
    context_ = context;
    this->adjust_formatter();
}

void formatter::PhysicalFmt::set_p_formatter(Fmt* p_formatter) {
    this->p_formatter_ = p_formatter;
    this->decorator_mode_ = true;
}

int formatter::Table::add_col(std::string new_title, std::vector<std::string> new_col) {
    this->ncol_++;
    this->titles_.push_back(new_title);
    this->data_.push_back(new_col);
    return this->ncol_-1;
}

void formatter::Table::adjust_col_width() {

    this->col_widths_.clear();
    int max_width = 0;
    for (int icol = 0; icol < this->ncol_; ++icol) {
        int s1 = this->titles_[icol].size();
        int s2 = 0;
        for (auto row : this->data_[icol]) {
            if (row.size() > s2) {
                s2 = row.size();
            }
        }
        if (s1 > s2) {
            if (flexible_width_) {
                this->col_widths_.push_back(s1);
            }
            else {
                if (s1 > max_width) {
                    max_width = s1;
                }
            }
        }
        else {
            if (flexible_width_) {
                this->col_widths_.push_back(s2);
            }
            else {
                if (s2 > max_width) {
                    max_width = s2;
                }
            }
        }
    }
    this->col_max_width_ = max_width;
    //initialize the vector col_widths_ by all the same number max_width
    for (int icol = 0; icol < this->ncol_; ++icol) {
        if (!flexible_width_) {
            this->col_widths_.push_back(max_width);
            this->total_width_ += max_width;
        }
        else {
            this->total_width_ += this->col_widths_[icol];
            if (this->col_max_width_ < this->col_widths_[icol]) {
                this->col_max_width_ = this->col_widths_[icol];
            }
        }
    }
    this->total_width_ += this->ncol_+1; // add the width of the delimiters
}

void formatter::Table::centerize_title() {
    for (int icol = 0; icol < this->ncol_; ++icol) {
        int s = this->titles_[icol].size();
        int n = this->col_widths_[icol];
        int n1 = (n-s)/2;
        int n2 = n-s-n1;
        std::string title = "";
        for (int i = 0; i < n1; ++i) {
            title += " ";
        }
        title += this->titles_[icol];
        for (int i = 0; i < n2; ++i) {
            title += " ";
        }
        this->titles_[icol] = title;
    }
}

void formatter::Table::reset() {
    this->overall_title = "";
    this->ncol_ = 0;

    this->col_delimiter_ = ' ';
    this->frame_switches_ = {1, 1, 0, 0};
    this->frame_delimiters_ = {'-', '-', '|', '|'};

    this->flexible_width_ = true;
    this->col_max_width_ = 0;
    this->total_width_ = 0;

    this->title_centered_ = false;
    this->clean();
}

void formatter::Table::clean() {
    this->titles_.clear();
    this->col_widths_.clear();
    this->data_.clear();
}

std::string formatter::Table::print_table() {

    this->adjust_col_width();
    if (this->title_centered_) {
        this->centerize_title();
    }
    std::stringstream ss;

    // the numbers of columns are allowed to be different, get the maximum number of rows
    int nrow_max = this->data_[0].size();
    std::vector<int> nrows;
    int ncol = this->ncol_;
    for (int icol = 0; icol < ncol; ++icol) {
        int s = this->data_[icol].size();
        nrows.push_back(s);
        if (nrow_max < s) {
            nrow_max = s;
        }
    }
    if (this->mode_ != 1) {
        if (this->overall_title.size() > 0) {
            ss << this->overall_title << std::endl;
        }
        // print the top frame
        if (this->frame_switches_[0]) {
            for (int iw = 0; iw < this->total_width_; ++iw) {
                ss << this->frame_delimiters_[0];
            }
            ss << std::endl;
        }
        // print the title
        if (this->frame_switches_[2]) {
            ss << this->frame_delimiters_[2];
        }
        for (int icol = 0; icol < ncol; ++icol) {
            ss << std::setw(this->col_widths_[icol]) 
            << std::setfill(' ') 
            << std::left 
            << this->titles_[icol];
            if (icol != ncol-1) {
                ss << this->col_delimiter_;
            }
        }
        if (this->frame_switches_[3]) {
            ss << this->frame_delimiters_[3];
        }
        ss << std::endl;
        // print the middle frame
        if (this->frame_mid_switch_) {
            for (int iw = 0; iw < this->total_width_; ++iw) {
                ss << this->frame_mid_delimiter_;
            }
            ss << std::endl;
        }
    }
    // print the data
    if (this->mode_ >= 0) {
        for (int irow = 0; irow < nrow_max; ++irow) {
            if (this->frame_switches_[2]) {
                ss << this->frame_delimiters_[2];
            }
            for (int icol = 0; icol < ncol; ++icol) {
                if (irow < nrows[icol]) {
                    ss << std::setw(this->col_widths_[icol]) 
                    << std::setfill(' ') 
                    << std::left 
                    << this->data_[icol][irow];
                }
                else {
                    ss << std::setw(this->col_widths_[icol]) 
                    << std::setfill(' ') 
                    << std::left 
                    << "";
                }
                if (icol != ncol-1) {
                    ss << this->col_delimiter_;
                }
            }
            if (this->frame_switches_[3]) {
                ss << this->frame_delimiters_[3];
            }
            ss << std::endl;
        }
    }
    // finally the bottom frame
    if (this->mode_ == 0) {
        if (this->frame_switches_[1]) {
            for (int iw = 0; iw < this->total_width_; ++iw) {
                ss << this->frame_delimiters_[1];
            }
            ss << std::endl;
        }
    }
    this->reset();
    return ss.str();
}

formatter::ContextFmt::ContextFmt() {
    this->p_phys_fmt_ = new formatter::PhysicalFmt("none", &(this->fmt_));
    // because the table is called inside ContextFmt, so the mode of the table is set to 'data'
    this->disable_title();
}

formatter::ContextFmt::~ContextFmt() {
    delete this->p_phys_fmt_;
}

void formatter::ContextFmt::set_context(std::string context) {
    if (strcmp(context.c_str(), this->context_.c_str()) != 0) {
        this->disable_iterative();
    }
    this->ncol_ = 0;
    this->context_ = context;
    this->phys_fmt_.clear();
    this->title_switch_ = 0;
    auto it = this->predefined_phys_fmt.find(context);
    if (it != this->predefined_phys_fmt.end()) {
        this->known_context_ = true;
        for (auto fmt : it->second.second) {
            this->phys_fmt_.push_back(fmt);
            this->ncol_++;
        }
        if (it->second.first == 1) {
            if (this->iterative_ >= 1) {
                this->enable_title();
            }
            else {
                this->disable_title();
            }
        }
        else if (it->second.first == 0) {
            // there is a title defined in this context, so it is needed to distinguish the title and the data
            this->enable_title();
        }
        else if (it->second.first == -1) {
            this->only_title();
        }
    }
    else {
        // do nothing
    }
}

void formatter::ContextFmt::set_context(int ncol, std::string* phys_fmt, int* nrows) {
    this->disable_iterative();
    std::string context = "none";
    this->context_ = context;
    this->ncol_ = ncol;
    this->nrows_.clear();
    this->phys_fmt_.clear();
    this->title_switch_ = 0;
    for (int icol = 0; icol < ncol; ++icol) {
        this->nrows_.push_back(nrows[icol]);
        this->phys_fmt_.push_back(phys_fmt[icol]);
    }
    auto it = this->predefined_phys_fmt.find(context);
    if (it != this->predefined_phys_fmt.end()) {
        this->known_context_ = true;
        if (it->second.first == 1) {
            this->disable_title();
        }
        else if (it->second.first == 0) {
            // there is a title defined in this context, so it is needed to distinguish the title and the data
            this->enable_title();
        }
        else if (it->second.first == -1) {
            this->only_title();
        }
    }
    else {
        // do nothing
    }
}

void formatter::ContextFmt::set_context(std::string context, int ncol, int* &nrows) {
    this->set_context(context);
    for (int i = 0; i < ncol; ++i) {
        this->nrows_.push_back(nrows[i]);
    }
}

void formatter::ContextFmt::set_context(std::vector<std::string> phys_fmt) {
    this->disable_iterative();
    std::string context = "none";
    this->context_ = context;
    this->ncol_ = phys_fmt.size();
    this->nrows_.clear();
    this->phys_fmt_.clear();
    this->title_switch_ = 0;
    for (int icol = 0; icol < this->ncol_; ++icol) {
        this->phys_fmt_.push_back(phys_fmt[icol]);
    }
    auto it = this->predefined_phys_fmt.find(context);
    if (it != this->predefined_phys_fmt.end()) {
        this->known_context_ = true;
        if (it->second.first == 1) {
            this->disable_title();
        }
        else if (it->second.first == 0) {
            // there is a title defined in this context, so it is needed to distinguish the title and the data
            this->enable_title();
        }
        else if (it->second.first == -1) {
            this->only_title();
        }
    }
    else {
        // do nothing
    }
}

template <typename T>
formatter::ContextFmt& formatter::ContextFmt::operator<<(const T& value) {

    if (this->title_switch_%2 == 0) {
        // it is a title
        this->cache_title_ = std::to_string(value);
    }
    else {
        // it is a data
        if (this->known_context_) {
            this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
        }
        else {
            this->p_phys_fmt_->set_context(this->default_phys_fmt_);
        }
        Table::add_col(this->cache_title_, (std::vector<std::string>){this->fmt_.format(value)});
        this->cache_title_ = "";
        this->icol_++;
    }
    if (Table::get_mode() == 1) {
        this->title_switch_ += 2;
    }
    else {
        this->title_switch_++;
    }
    return *this;
}
template <>
formatter::ContextFmt& formatter::ContextFmt::operator<<(const std::string& value) {

    if (this->title_switch_%2 == 0) {
        // it is a title
        this->cache_title_ = value;  
    }
    else {
        // it is a data
        if (this->known_context_) {
            // known context, will be formatted
            this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
        }
        else {
            // unknown context, will not be formatted
            this->p_phys_fmt_->set_context("none");
        }
        Table::add_col(this->cache_title_, (std::vector<std::string>){this->fmt_.format(value)});
        this->icol_++;
    }
    if (Table::get_mode() == 1) {
        this->title_switch_ += 2;
    }
    else {
        this->title_switch_++;
    }
    return *this;
}

formatter::ContextFmt& formatter::ContextFmt::operator<< (char const* value)
{
    std::string value_(value);
    return *this << value_;
}
template formatter::ContextFmt& formatter::ContextFmt::operator<< <double>(double const& value);
template formatter::ContextFmt& formatter::ContextFmt::operator<< <int>(int const& value);

template <typename T>
formatter::ContextFmt& formatter::ContextFmt::operator<<(const std::vector<T>& value) {
    // well if a whole vector is input, it is definitely a data... or presently I can't think of a case that a vector is a title
    std::vector<std::string> value_;
    if (this->known_context_) {
        this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
    }
    else if (this->icol_ < this->phys_fmt_.size()) {
        this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
    }
    else {
        this->p_phys_fmt_->set_context(this->default_phys_fmt_);
    }

    for (auto v : value) {
        value_.push_back(this->fmt_.format(v));
    }
    Table::add_col(this->cache_title_, value_);
    this->cache_title_ = "";
    this->icol_++;
    this->title_switch_ += 2;
    return *this;
}
template formatter::ContextFmt& formatter::ContextFmt::operator<< <double>(std::vector<double> const& value);
template formatter::ContextFmt& formatter::ContextFmt::operator<< <std::string>(std::vector<std::string> const& value);
template formatter::ContextFmt& formatter::ContextFmt::operator<< <int>(std::vector<int> const& value);
template formatter::ContextFmt& formatter::ContextFmt::operator<< <char>(char const& value);

template <typename T>
formatter::ContextFmt& formatter::ContextFmt::operator<<(T*& value) {
    std::vector<std::string> value_;
    if (this->known_context_) {
        this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
    }
    else if (this->icol_ < this->phys_fmt_.size()) {
        this->p_phys_fmt_->set_context(this->phys_fmt_[this->icol_]);
    }
    else {
        this->p_phys_fmt_->set_context(this->default_phys_fmt_);
    }
    for (int i = 0; i < this->nrows_[this->icol_]; ++i) {
        value_.push_back(this->fmt_.format(value[i]));
    }
    Table::add_col(this->cache_title_, value_);
    this->cache_title_ = "";
    this->icol_++;
    this->title_switch_ += 2;
    return *this;
}
template formatter::ContextFmt& formatter::ContextFmt::operator<< <double>(double*& value);
template formatter::ContextFmt& formatter::ContextFmt::operator<< <int>(int*& value);
template formatter::ContextFmt& formatter::ContextFmt::operator<< <std::string>(std::string*& value);
template formatter::ContextFmt& formatter::ContextFmt::operator<< <char>(char*& value);

void formatter::ContextFmt::reset() {
    Table::reset();
    this->phys_fmt_.clear();
    this->cache_title_ = "";
    this->icol_ = 0;
    this->ncol_ = 0;
    this->nrows_.clear();
    this->title_switch_ = 0;
    this->set_context(this->context_);
}

std::string formatter::ContextFmt::str(bool next_line) {
    if (this->iterative_ == 1)
    {
        Table::set_mode(0);
        Table::disable_down_frame();
        this->iterative_ += 1;
    }
    else if (this->iterative_ > 1)
    {
        Table::set_mode(1);
        this->iterative_ += 1;
    }
    std::string str = this->print_table();
    if (!next_line) {
        str.pop_back();
    }
    this->reset();
    return str;
}

void formatter::ContextFmt::print_status() const {
    std::cout << "-----------------status-----------------" << std::endl;
    std::cout << "context: " << this->context_ << std::endl;
    std::cout << "known_context: " << this->known_context_ << std::endl;
    std::cout << "ncol: " << this->ncol_ << std::endl;
    std::cout << "nrows: ";
    for (auto n : this->nrows_) {
        std::cout << n << " ";
    }
    std::cout << std::endl;
    std::cout << "phys_fmt: ";
    for (auto fmt : this->phys_fmt_) {
        std::cout << fmt << " ";
    }
    std::cout << std::endl;
    std::cout << "cache_title: " << this->cache_title_ << std::endl;
    std::cout << "icol: " << this->icol_ << std::endl;
    std::cout << "title_switch: " << this->title_switch_ << std::endl;
    std::cout << "iterative: " << this->iterative_ << std::endl;
    std::cout << "fmt: " << std::endl;
    std::cout << "width: " << this->fmt_.get_width() << std::endl;
    std::cout << "precision: " << this->fmt_.get_precision() << std::endl;
    std::cout << "fillChar: " << this->fmt_.get_fillChar() << std::endl;
    std::cout << "fixed: " << this->fmt_.get_fixed() << std::endl;
    std::cout << "right: " << this->fmt_.get_right() << std::endl;
    std::cout << "error: " << this->fmt_.get_error() << std::endl;
    std::cout << "p_phys_fmt: " << this->p_phys_fmt_ << std::endl;
    std::cout << "mode: " << Table::get_mode() << std::endl;
    std::cout << "col_delimiter: " << Table::get_col_delimiter() << std::endl;
    std::cout << "frame_switches: ";
    for (auto s : Table::get_frame_switches()) {
        std::cout << s << " ";
    }
    std::cout << std::endl;
    std::cout << "frame_delimiters: ";
    for (auto d : Table::get_frame_delimiters()) {
        std::cout << d << " ";
    }
    std::cout << std::endl;
    std::cout << "flexible_width: " << Table::get_flexible_width() << std::endl;
    std::cout << "frame_mid_switch: " << Table::get_frame_mid_switch() << std::endl;
    std::cout << "frame_mid_delimiter: " << Table::get_frame_mid_delimiter() << std::endl;
    std::cout << "col_widths: ";
    for (auto w : Table::get_col_widths()) {
        std::cout << w << " ";
    }
    std::cout << std::endl;
    std::cout << "col_max_width: " << Table::get_col_max_width() << std::endl;
    std::cout << "total_width: " << Table::get_total_width() << std::endl;
    std::cout << "----------------------------------------" << std::endl;
}

formatter::ContextFmt context;
