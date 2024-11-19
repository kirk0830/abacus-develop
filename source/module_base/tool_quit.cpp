#include "module_base/tool_quit.h"

#ifdef __MPI
#include "mpi.h"
#include "module_base/parallel_global.h"
#include "module_base/parallel_comm.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __NORMAL
#else
#include "module_base/global_variable.h"
#include "module_parameter/parameter.h"
#include "module_base/global_file.h"
#include "module_base/timer.h"
#include "module_base/memory.h"
#endif

namespace ModuleBase
{
//==========================================================
// GLOBAL FUNCTION :
// NAME : WARNING( write information into GlobalV::ofs_warning)
// NAME : QUIT( exit the running program)
// NAME : WARNING_QUIT( write information into
//           GlobalV::ofs_warning , and then quit)
//==========================================================
void WARNING(const std::string &file,const std::string &description)
{
#ifdef __NORMAL /* what is this??? */
#else
    if (GlobalV::MY_RANK==0)
    {
        GlobalV::ofs_warning << " " << file <<"  warning : "<< description<<std::endl;
    }
#endif
    return;
}

void QUIT()
{
	/* it is quite strange to exit with code 0 by default... */
    QUIT(0);
}

void QUIT(const int exitcode)
{
#ifdef __NORMAL /* what is this??? */
#else
    ModuleBase::timer::finish(GlobalV::ofs_running , !GlobalV::MY_RANK);
    ModuleBase::Global_File::close_all_log(GlobalV::MY_RANK);
    std::cout<<" See output information in : "<<PARAM.globalv.global_out_dir<<std::endl;
#endif
#ifdef _OPENMP // merge all threads of one process into one thread
    if (omp_in_parallel())
    {
        omp_set_num_threads(1);
        std::cout << "Terminating ABACUS with multithreading environment." << std::endl;
    }
    assert(!omp_in_parallel()); /* avoid the case that death thread calls fork() */
#endif
#ifdef __MPI /* if it is MPI run, finalize first, then exit */
    std::cout << "Terminating ABACUS with multiprocessing environment." << std::endl;
    Parallel_Global::finalize_mpi();
    /* but seems this is the only correct way to terminate the MPI */
#endif
    exit(exitcode);
}

void WARNING_QUIT(const std::string &file, const std::string &description)
{
    WARNING_QUIT(file, description, 0); 
	/* really? we return with 0? it is something like "we exit the program normally" */
}

void WARNING_QUIT(const std::string &file, const std::string &description, int exitcode)
{
    const std::string banner = ""
    " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
    "                         NOTICE                           \n"
    " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
#ifdef __NORMAL /* what is this??? */
    std::cout << banner << std::endl;
#else
    std::cout << banner.size() << std::endl;
    const int max_width = banner.size() / 3 - 1; // minus 1 because of the "\n"
    // wrap the description
    std::string wrapped_desc = "";
    std::string::size_type pos = 0;
    while (pos < description.size())
    {
        wrapped_desc += " ";
        wrapped_desc += description.substr(pos, max_width - 2);
        // because the leading whitespace and tailing "\n"
        wrapped_desc += "\n";
        pos += max_width - 2;
    }
    const std::string warnmsg = "\n"
    " " + description + "\n"
    " CHECK IN FILE : " + PARAM.globalv.global_out_dir + "warning.log\n\n";
    std::cout << "\n" << banner << warnmsg << banner << std::endl;
    GlobalV::ofs_running << "\n" << banner << warnmsg << banner << std::endl;
    WARNING(file,description);
    GlobalV::ofs_running<<" Check in file : "<<PARAM.globalv.global_out_dir<<"warning.log"<<std::endl;
#endif
    QUIT(exitcode);
}

//Check and print warning information for all cores.
//Maybe in the future warning.log should be replaced by error.log.
void CHECK_WARNING_QUIT(const bool error_in, const std::string &file,const std::string &calculation,const std::string &description)
{
#ifdef __NORMAL
// only for UT, do nothing here
#else
    if(error_in)
    {
        //All cores will print inforamtion
        std::cout.clear();
        if(!GlobalV::ofs_running.is_open())
        {
            std::string logfile = PARAM.globalv.global_out_dir + "running_" + calculation + ".log";
            GlobalV::ofs_running.open( logfile.c_str(), std::ios::app );
        }
        if(!GlobalV::ofs_warning.is_open())
        {
            std::string warningfile = PARAM.globalv.global_out_dir + "warning.log";
            GlobalV::ofs_warning.open( warningfile.c_str(), std::ios::app );
        }

        //print error information
        std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout << " ERROR! " << description << std::endl;
        std::cout << " CHECK IN FILE : " << PARAM.globalv.global_out_dir << "warning.log" << std::endl;
        std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        GlobalV::ofs_running << " ERROR! CHECK IN FILE : " << PARAM.globalv.global_out_dir << "warning.log" << std::endl;
        GlobalV::ofs_warning << std::endl;
        GlobalV::ofs_warning << " ERROR! " << file << ", core " << GlobalV::MY_RANK+1 << ": " << description << std::endl;
        GlobalV::ofs_warning << std::endl;
        exit(1);
    }
#endif
    return;
}

}
