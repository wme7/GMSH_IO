#ifndef PARAMETERS_CLASS
#define PARAMETERS_CLASS

#include "Globals.hpp"
#include "INIreader.hpp"

// Read parameters form an standard INIfile

Parameters INIparser(std::string ini_file)
{
    // Initialize Parameters object
    Parameters m;

    // Read the ini_file
    INIReader reader(ini_file);

    // If we can't load the ini_file
    if (reader.ParseError() < 0) {
        std::cout << "Can't load " << ini_file << std::endl;
        std::exit(-1);
    }

    // Parse the parameters for simulation as:
    // m.variable = reader.Get("section","field","default_value");

    // [run]
    m.m_solver  = reader.Get("run","solver","none");
    m.m_mesher  = reader.Get("run","mesher","none");
    m.m_binary  = reader.Get("run","binary","none");
    m.m_name    = reader.Get("run","name","none");
    if (m.m_name == "none")
    {
        size_t lastindex = ini_file.find_last_of("."); 
        std::string rawname = ini_file.substr(0, lastindex);
        m.m_name = rawname;
    }
    
    //[mesh]
    if ((m.m_mesher != "none") | (m.m_binary != "none") | (m.m_solver != "none"))
    {
        m.m_N         = reader.GetSize_t("mesh","N",1);
        m.m_mesh_name = reader.Get("mesh","name","none");
        if (m.m_mesh_name == "none")
        {
            std::cout << "ERROR: no mesh file has been provided!" << std::endl; 
            std::exit(-1);
        }
    }
    
    //[time]
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_t_start = reader.GetReal("time","start",0.0);
        m.m_dt_user = reader.GetReal("time","dt",0.0);
        m.m_Nt_user = reader.GetSize_t("time","Nt",0);
        m.m_t_end = reader.GetReal("time","end",m.m_dt_user * m.m_Nt_user); // follow the user request: t_end = dt*Nt.
        if (m.m_t_end == 0.) {
            std::cout << "WARNING: output time not set directly (" << m.m_t_start << " < t < t_end? )" << std::endl;
        } 
        m.m_sv = reader.GetSize_t("time","save",10);
        m.m_CFL = reader.GetReal("time","CFL",0.5);
        if ( (m.m_CFL <= 0.0) | (m.m_CFL > 1.0) )
        {
            std::cout << "WARNING: CFL condition should be CFL < 1. Setting default value CFL = 0.5" << std::endl;
            m.m_CFL = 0.5;
        }
    }

    //[initial]
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_initial_mode = reader.Get("initial","mode","homogeneous");
        std::cout<< "IC ......... : " << m.m_initial_mode << std::endl;

        if (m.m_initial_mode == "homogeneous")
        {
            m.m_x0 = 0.;
            m.m_y0 = 0.;
            m.m_z0 = 0.;
        }
        if ((m.m_initial_mode == "gaussian") | (m.m_initial_mode == "gen"))
        {
            m.m_Amp = reader.GetReal("initial","A",10.0);
            m.m_x0  = reader.GetReal("initial","x0",0.0);
            m.m_y0  = reader.GetReal("initial","y0",0.0);
            m.m_z0  = reader.GetReal("initial","z0",0.0);
            m.m_Sgm = reader.GetReal("initial","sigma",0.3*1500/1E5);
        }
        if (m.m_initial_mode == "plane") 
        {
            m.m_Amp   = reader.GetReal("initial","A",10.0);
            m.m_theta = reader.GetReal("initial","theta",10.0);
            m.m_R     = reader.GetReal("initial","R",0.02);
            m.m_x0    = reader.GetReal("initial","x0",0.0);
            m.m_y0    = reader.GetReal("initial","y0",0.0);
            m.m_z0    = reader.GetReal("initial","z0",0.0);
            m.m_Sgm   = reader.GetReal("initial","sigma",0.02);
        }
        if (m.m_initial_mode == "read") 
        {
            // do nothing !
        }
    }
    
    //[boundaries]
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_boundaries_mode = reader.Get("boundaries","mode","auto");
        std::cout<< "BC ......... : " << m.m_boundaries_mode << std::endl;
        
        if(m.m_boundaries_mode=="auto")
        {
            // Do nothing and just follow the BCs set in the meshfile.
            std::cout<< "BC ......... : " << "get from mesh"  << std::endl;
        }
        if (m.m_boundaries_mode=="gaussian_spot")
        {
            std::cout<< "BC gaussian spot detected : "  << std::endl;
            m.m_BCSSfile              = reader.Get("boundaries","file","BC_signals.h5");
            m.m_nb_boundaries_sources = reader.GetSize_t("boundaries","nb",0);
        }
        if (m.m_boundaries_mode=="impose")
        {
            // Follow the BCs set in the meshfile. However, 
            // in regions with impossed (pressure and/or velocity) conditions, we set:
            m.m_p0_BC = reader.GetReal("boundaries","p0",0.0);
            m.m_u0_BC = reader.GetReal("boundaries","u0",0.0);
            m.m_v0_BC = reader.GetReal("boundaries","v0",0.0);
            m.m_w0_BC = reader.GetReal("boundaries","w0",0.0);
            m.m_f0_BC = reader.GetReal("boundaries","f0",0.0);
        }
    }

    //[fields] : Query selected solution fields throught out the computation
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_query_fields_mode = reader.Get("fields","mode","none");
        std::cout<< "QQ ......... : " << m.m_query_fields_mode << std::endl;

        if (m.m_query_fields_mode=="true")
        {
            m.m_tq_start = reader.GetReal("fields","start",m.m_t_start); // Set equal to t_start
            m.m_tq_end = reader.GetReal("fields","end",m.m_t_end); // Set equal to t_end
            
            m.m_qmax_eq = reader.GetInteger("fields","qmax_eq",-1);
            m.m_qmax_0  = reader.GetReal   ("fields","qmax_0",0.0);
            m.m_qmin_eq = reader.GetInteger("fields","qmin_eq",-1);
            m.m_qmin_0  = reader.GetReal   ("fields","qmin_0",0.0);
            m.m_qsum_eq = reader.GetInteger("fields","qsum_eq",-1);
            m.m_qsum_0  = reader.GetReal   ("fields","qsum_0",0.0);
            m.m_qsum2_eq= reader.GetInteger("fields","qsum2_eq",-1);
            m.m_qsum2_0 = reader.GetReal   ("fields","qsum2_0",0.0);
            m.m_qsum3_eq= reader.GetInteger("fields","qsum3_eq",-1);
            m.m_qsum3_0 = reader.GetReal   ("fields","qsum3_0",0.0);
            m.m_qsum4_eq= reader.GetInteger("fields","qsum4_eq",-1);
            m.m_qsum4_0 = reader.GetReal   ("fields","qsum4_0",0.0);
            m.m_qSkewness_eq = reader.GetInteger("fields","qSkewness_eq",-1);
            m.m_qSkewness_0  = reader.GetReal   ("fields","qSkewness_0",0.0);
            m.m_qKurtosis_eq = reader.GetInteger("fields","qKurtosis_eq",-1);
            m.m_qKurtosis_0  = reader.GetReal   ("fields","qKurtosis_0",0.0);
            m.m_qIntensity_eq= reader.GetInteger("fields","qIntensity_eq",-1);
            m.m_qIntensity_0 = reader.GetReal   ("fields","qIntensity_0",0.0);
            m.m_qVelocity_eq = reader.GetInteger("fields","qVelocity_eq",-1);
            m.m_qVelocity_0  = reader.GetReal   ("fields","qVelocity_0",0.0);
            m.m_qRMS_eq = reader.GetInteger("fields","qRMS_eq",-1);
            m.m_qRMS_0  = reader.GetReal   ("fields","qRMS_0",0.0);
        } else {
            m.m_qmax_eq =-1; m.m_qmax_0 = 0.;
            m.m_qmin_eq =-1; m.m_qmin_0 = 0.;
            m.m_qsum_eq =-1; m.m_qsum_0 = 0.;
            m.m_qsum2_eq =-1; m.m_qsum2_0 = 0.;
            m.m_qsum3_eq =-1; m.m_qsum3_0 = 0.;
            m.m_qsum4_eq =-1; m.m_qsum4_0 = 0.;
            m.m_qSkewness_eq =-1; m.m_qSkewness_0 = 0.;
            m.m_qKurtosis_eq =-1; m.m_qKurtosis_0 = 0.;
            m.m_qVelocity_eq =-1; m.m_qVelocity_0 = 0.;
            m.m_qIntensity_eq =-1; m.m_qIntensity_0 = 0.;
            m.m_qRMS_eq =-1; m.m_qRMS_0 = 0.;
        }
    }

    //[parameters]: Physical Parameters
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_physical_parameter_mode = reader.Get("parameters","mode","default");
        std::cout<< "PP ......... : " << m.m_physical_parameter_mode << " : ";

        if (m.m_physical_parameter_mode=="default")
        {
            if((m.m_solver=="ParadigmS2D_acouSolver.run") | (m.m_solver=="ParadigmS3D_acouSolver.run"))
            {
                m.m_Npara = 2;
                m.m_physical_parameter_names.push_back("c0");
                m.m_physical_parameter_values.push_back(1500.);
                m.m_physical_parameter_names.push_back("rho0");
                m.m_physical_parameter_values.push_back(1000.);
                std::cout << "homogeneous water" << std::endl;
            }
            if((m.m_solver=="ParadigmS2D_elastoSolver.run") | (m.m_solver=="ParadigmS3D_elastoSolver.run"))
            {
                m.m_Npara = 3;
                m.m_physical_parameter_names.push_back("lambda0");
                m.m_physical_parameter_values.push_back(60080000000.);
                m.m_physical_parameter_names.push_back("mu0");
                m.m_physical_parameter_values.push_back(20325000000.);
                m.m_physical_parameter_names.push_back("rho0");
                m.m_physical_parameter_values.push_back(2700.);
                std::cout << "homogeneous solid" << std::endl;
            }
            if(m.m_solver=="ParadigmS2D_LEE_Solver.run")
            {
                m.m_Npara = 5;
                m.m_physical_parameter_names.push_back("c0");
                m.m_physical_parameter_values.push_back(340.);
                m.m_physical_parameter_names.push_back("rho0");
                m.m_physical_parameter_values.push_back(1.2);
                m.m_physical_parameter_names.push_back("u0");
                m.m_physical_parameter_values.push_back(0.0);
                m.m_physical_parameter_names.push_back("v0");
                m.m_physical_parameter_values.push_back(0.0);
                m.m_physical_parameter_names.push_back("p0");
                m.m_physical_parameter_values.push_back(100000.);
                std::cout << "homogeneous standard air" << std::endl;
            }
            if (m.m_solver=="ParadigmS3D_LEE_Solver.run")
            {
                m.m_Npara = 6;
                m.m_physical_parameter_names.push_back("c0");
                m.m_physical_parameter_values.push_back(340.);
                m.m_physical_parameter_names.push_back("rho0");
                m.m_physical_parameter_values.push_back(1.2);
                m.m_physical_parameter_names.push_back("u0");
                m.m_physical_parameter_values.push_back(0.0);
                m.m_physical_parameter_names.push_back("v0");
                m.m_physical_parameter_values.push_back(0.0);
                m.m_physical_parameter_names.push_back("w0");
                m.m_physical_parameter_values.push_back(0.0);
                m.m_physical_parameter_names.push_back("p0");
                m.m_physical_parameter_values.push_back(100000.);
                std::cout << "homogeneous standard air" << std::endl;
            }
            if (m.m_solver=="none") 
            {
                std::cout << "no solver has been set!" << std::endl;
            }
        }
        if ((m.m_physical_parameter_mode=="gen") | (m.m_physical_parameter_mode=="homogeneous"))
        {
            // Get number of parameters
            m.m_Npara = reader.GetSize_t("parameters","Npara",0);

            for(size_t n=1; n<=m.m_Npara; n++)
            {           
                // Get Parameters names
                std::string name = reader.Get("parameters","name_para"+std::to_string(n),"none");
                m.m_physical_parameter_names.push_back(name);

                // Get Parameters values
                double value = reader.GetReal("parameters","value_para"+std::to_string(n),0.0);
                m.m_physical_parameter_values.push_back(value);
            }
            std::cout << std::endl;
        }
        if (m.m_physical_parameter_mode=="read") // normaly for domains with inhomogeneous physical properties
        {
            std::cout << " from file." << std::endl;

            // Number of parameters
            m.m_Npara = 0;
        }
    }
    
    //[sources]
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_sources_mode = reader.Get("sources","mode","false");
        m.m_nb_sources   = reader.GetSize_t("sources","nb",0);
        std::cout<< "SS ......... : " << m.m_sources_mode << " : INIfile has " << m.m_nb_sources << std::endl;

        if ((m.m_sources_mode == "gen") | (m.m_sources_mode == "gaussian") | (m.m_sources_mode == "true"))
        {
            for(size_t ns=1; ns<=m.m_nb_sources; ns++)
            {
                m.m_A.push_back(        reader.GetReal("sources",  "A"   + std::to_string(ns),0) );                
                m.m_sigma.push_back(    reader.GetReal("sources","sigma" + std::to_string(ns),0) );
                m.m_xs.push_back(       reader.GetReal("sources",  "xs"  + std::to_string(ns),0) );
                m.m_ys.push_back(       reader.GetReal("sources",  "ys"  + std::to_string(ns),0) );
                m.m_zs.push_back(       reader.GetReal("sources",  "zs"  + std::to_string(ns),0) );
                m.m_f0.push_back(       reader.GetReal("sources",  "f0"  + std::to_string(ns),0) );
                m.m_phi.push_back(      reader.GetReal("sources",  "phi" + std::to_string(ns),0) );
                m.m_T0.push_back(       reader.GetReal("sources",  "T0"  + std::to_string(ns),0) );
                m.m_sources_eq.push_back(reader.GetSize_t("sources","eq" + std::to_string(ns),0) );
            }
        }
        if (m.m_sources_mode == "read")
        {
            m.m_SSfile = reader.Get("sources","file","Signals.h5");
            for(size_t ns=1; ns<=m.m_nb_sources; ns++)
            {
                m.m_A.push_back(        reader.GetReal("sources",  "A"   + std::to_string(ns),0) );                
                m.m_sigma.push_back(    reader.GetReal("sources","sigma" + std::to_string(ns),0) );
                m.m_xs.push_back(       reader.GetReal("sources",  "xs"  + std::to_string(ns),0) );
                m.m_ys.push_back(       reader.GetReal("sources",  "ys"  + std::to_string(ns),0) );
                m.m_zs.push_back(       reader.GetReal("sources",  "zs"  + std::to_string(ns),0) );
                m.m_sources_eq.push_back(reader.GetSize_t("sources","eq" + std::to_string(ns),0) );
            }
        }
        if ((m.m_sources_mode == "multiread") | (m.m_sources_mode == "gaussian_multiread"))
        {
            m.m_SSfile = reader.Get("sources","file","Signals.h5");
            std::cout<< "NOTE: the data and parameters for " << m.m_nb_sources << " source(s) will be obtain from " << m.m_SSfile << std::endl;
        }
    }

    //[output]
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_xmf_mode = reader.Get("output","xmf","true"); // true: output xmf files
        m.m_ipynb_mode = reader.Get("output","ipynb","true"); // true: output ipynb solution
        m.m_terminal_mode = reader.Get("output","terminal","true"); // true: print all terminal output
    }

    //[record]
    if ((m.m_solver != "none") | (m.m_binary != "none"))
    {
        m.m_rec_mode = reader.Get("record","mode","all_fields"); // record by default "all_fields"
        if (m.m_rec_mode == "field")
        {
            m.m_rec_eq = reader.GetSize_t("record","eq",0); // set and specific field/equation to be recorded!
        } 
        if ((m.m_rec_mode != "field") & (m.m_rec_mode != "all_fields"))
        {
            m.m_rec_mode = "off"; // Cancel rec-objects in the mesh. Do not record anything!
        }
    }

    // If everything is ok
    return m;
}

#endif
