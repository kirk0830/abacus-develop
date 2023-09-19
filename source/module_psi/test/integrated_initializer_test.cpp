#include <iostream>
#include <fstream>
#include <string>
#include <regex>

#include <gtest/gtest.h>

class IntegratedInitializerTest : public ::testing::Test {
protected:
    int error = 0;
    std::ifstream ifs_new;
    std::ifstream ifs_old;
    int n_match_max;
    int n_match;

    std::regex pattern;
    virtual void SetUp() {
        std::cout << "IntegratedInitializerTest SetUp" << std::endl;
        this->n_match = 0;
        this->n_match_max = 10;
        this->pattern = std::regex("\\((-?[0-9]+\\.[0-9]*),(-?[0-9]+\\.[0-9]*)\\)");
        this->error = std::system("cp ../../../../source/module_psi/test/support/KPT ./KPT");
        this->error = std::system("cp ../../../../source/module_psi/test/support/STRU ./STRU");
        this->error = std::system("cp ../../../../source/module_psi/test/support/Si_NCSR_ONCVPSP_v0.5_dojo.upf ./Si_NCSR_ONCVPSP_v0.5_dojo.upf");
        this->error = std::system("cp ../../../../source/module_psi/test/support/Si_gga_8au_60Ry_2s2p1d.orb ./Si_gga_8au_60Ry_2s2p1d.orb");
    }

    virtual void TearDown() {
        std::cout << "IntegratedInitializerTest TearDown" << std::endl;
        this->error = std::system("rm ./psig_0_kpt.out");
        // finalize
        this->error = std::system("rm ./Si_NCSR_ONCVPSP_v0.5_dojo.upf");
        this->error = std::system("rm ./Si_gga_8au_60Ry_2s2p1d.orb");
        this->error = std::system("rm ./KPT");
        this->error = std::system("rm ./STRU");
        this->error = std::system("rm -rf ./OUT*");
        this->error = std::system("rm ./time.json");
    }
};

TEST_F(IntegratedInitializerTest, CalPsiGRandom) {

    int error = 0;
    // init_wfc = "random"
    error = std::system("cp ../../../../source/module_psi/test/support/random_new ./INPUT");
    error = std::system("../../../abacus");
    error = std::system("rm ./INPUT");
    error = std::system("cp ../../../../source/module_psi/test/support/random_old ./INPUT");
    error = std::system("../../../abacus");
    error = std::system("rm ./INPUT");
    // compare wavefunctions
    this->ifs_new.open("./psig_0_kpt.out");
    this->ifs_old.open("./psig_0_kpt_old.out");
    std::string line_new, line_old;
    std::smatch match_new, match_old;

    while (ifs_new >> line_new && ifs_old >> line_old) {
        if (std::regex_search(line_new, match_new, this->pattern) && std::regex_search(line_old, match_old, pattern)) {
            EXPECT_EQ(match_new.size(), match_old.size());
            double real_new = std::stod(match_new[1]);
            double imag_new = std::stod(match_new[2]);
            double real_old = std::stod(match_old[1]);
            double imag_old = std::stod(match_old[2]);
            EXPECT_NEAR(real_new, real_old, 1e-10);
            EXPECT_NEAR(imag_new, imag_old, 1e-10);
            this->n_match++;
            if (this->n_match == this->n_match_max) break;
        }
    }
    this->ifs_new.close();
    this->ifs_old.close();

    error = std::system("rm ./psig_0_kpt_old.out");
}

TEST_F(IntegratedInitializerTest, CalPsiGAtomic) {

    int error = 0;
    std::string line_new;
    std::smatch match_new;
    // compare wavefunctions
    // init_wfc = "atomic"
    error = std::system("cp ../../../../source/module_psi/test/support/atomic_new ./INPUT");
    error = std::system("../../../abacus");
    error = std::system("rm ./INPUT");
    // compare wavefunctions with reference value
    std::vector<double> reals_ref = {
        -0.0000060004, 0.0000170688, 0.0000537252, -0.0006497418, -0.0050094394,
        -0.0129273778, 0.1125669791, 0.9193169677, 0.1125669791, -0.0129273778
    };
    std::vector<double> imags_ref = {
        0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -0.0, -0.0, -0.0
    };
    double* reals_read = new double[reals_ref.size()];
    double* imags_read = new double[imags_ref.size()];
    this->ifs_new.open("./psig_0_kpt.out");

    bool find_psi = false;
    while (ifs_new >> line_new) {
        if (line_new == "atomic") {
            std::cout << "Have find correct psi" << std::endl;
            find_psi = true;
            break;
        };
    }
    ASSERT_TRUE(find_psi);
    while (ifs_new >> line_new) {
        if (regex_search(line_new, match_new, this->pattern)) { //check if the line contains a complex number tuple
            double real = stod(match_new[1]); //extract the real part of the complex number
            double imag = stod(match_new[2]); //extract the imaginary part of the complex number
            reals_read[this->n_match] = real;
            imags_read[this->n_match] = imag;
            this->n_match++;
            if (this->n_match == this->n_match_max) break;
        }
    }
    this->ifs_new.close();
    for (int i = 0; i < reals_ref.size(); i++) {
        EXPECT_NEAR(reals_read[i], reals_ref[i], 1e-5);
        EXPECT_NEAR(imags_read[i], imags_ref[i], 1e-5);
    }
    delete[] reals_read;
    delete[] imags_read;
}

TEST_F(IntegratedInitializerTest, CalPsiGNao) {

    int error = 0;
    std::string line_new;
    std::smatch match_new;
    // compare wavefunctions
    // init_wfc = "nao"
    error = std::system("cp ../../../../source/module_psi/test/support/nao_new ./INPUT");
    error = std::system("../../../abacus");
    error = std::system("rm ./INPUT");
    std::vector<double> reals_ref = {
        0.0000005904, -0.0000000000, 0.0000083311, -0.0005825703, -0.0039503908,
        0.0000000000, 0.0477907663, 0.8852328311, 0.1021402963, 0.0000000000
    };
    std::vector<double> imags_ref = {
        -0.0000002762, -0.0000000000, 0.0000178056, -0.0002111799, 0.0018483616,
        -0.0000000000, 0.1021402963, 0.3208940657, -0.0477907663, -0.0000000000
    };
    double* reals_read = new double[reals_ref.size()];
    double* imags_read = new double[imags_ref.size()];
    // compare wavefunctions with reference value
    this->ifs_new.open("./psig_0_kpt.out");

    bool find_psi = false;
    while (ifs_new >> line_new) {
        if (line_new == "nao") {
            std::cout << "Have find correct psi" << std::endl;
            find_psi = true;
            break;
        };
    }
    ASSERT_TRUE(find_psi);

    while (ifs_new >> line_new) {
        if (regex_search(line_new, match_new, this->pattern)) { //check if the line contains a complex number tuple
            double real = stod(match_new[1]); //extract the real part of the complex number
            double imag = stod(match_new[2]); //extract the imaginary part of the complex number
            reals_read[this->n_match] = real;
            imags_read[this->n_match] = imag;
            n_match++;
            if (this->n_match == this->n_match_max) break;
        }
    }
    this->ifs_new.close();
    for (int i = 0; i < reals_ref.size(); i++) {
        EXPECT_NEAR(reals_read[i], reals_ref[i], 1e-5);
        EXPECT_NEAR(imags_read[i], imags_ref[i], 1e-5);
    }
    delete[] reals_read;
    delete[] imags_read;
}

#include <unistd.h>
int main(int argc, char **argv)
{
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    std::cout << "present directory: " << cwd << std::endl;
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
	return result;
}