// Extension of the program that manages the input and output.

#include "io.h"

Settings IO::loadSettings(const char *configFile) {
  Config cfg;
  Settings settings;

  // Read the file. If there is an error, report it and exit -------------------
  try
  {
    cfg.readFile(configFile);
  }
  catch(const FileIOException &fioex)
  {
    cerr << "I/O error while reading file." << std::endl;
  }
  catch(const ParseException &pex)
  {
    cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
         << " - " << pex.getError() << endl;
  }

  // Get the Project name ------------------------------------------------------
  try
  {
    string name = cfg.lookup("name");
    cout << "Project name: " << name << endl << endl;
    settings.name = name;
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No 'name' setting in configuration file." << endl;
  }

  const Setting &root = cfg.getRoot();

  // Get the Data Set information ----------------------------------------------
  try
  {
    const Setting &dataSet = root["dataSet"];

    // Data Set type
    try
    {
      string type = dataSet.lookup("type");
      settings.type = type;
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'dataSet type' setting in 'dataSet' section in configuration file." << endl;
    }

    // Data Set generation
    if (settings.type == "generate") {
      dataSet.lookupValue("name", settings.testCase.name);

      if (dataSet.lookupValue("model", settings.testCase.model)) {
        dataSet.lookupValue("nParam", settings.testCase.nParam);

        const Setting &muList = root["dataSet"]["mu"];
        int muParams = muList.getLength();

        for (int i = 0; i < muParams; i++) {
          int muModes = muList[i].getLength();
          vector<double> muParam (muModes);

          for (int j = 0; j < muModes; j++) {
            double mu = muList[i][j];
            muParam[j] = mu;
          }
          settings.testCase.mu.push_back(muParam);
        }

        const Setting &sigmaList = root["dataSet"]["sigma"];
        int sigmaParams = sigmaList.getLength();

        for (int i = 0; i < sigmaParams; i++) {
          int sigmaModes = sigmaList[i].getLength();
          vector<double> sigmaParam (sigmaModes);

          for (int j = 0; j < sigmaModes; j++) {
            double sigma = sigmaList[i][j];
            sigmaParam[j] = sigma;
          }
          settings.testCase.sigma.push_back(sigmaParam);
        }
      }

      // Get Charactieristic Diameter vector
      try
      {
        const Setting &D = root["dataSet"]["D"];

        bool generate;
        if (D.lookupValue("generate", generate)
            && generate == true) {

          double first = D["vec"][0];
          double step  = D["vec"][1];
          double last  = D["vec"][2];

          for (double value = first; value <= last; value += step) {
            settings.testCase.D.push_back(value);
          }
        }
        else
        {
          int count = D["vec"].getLength();
          settings.testCase.D.resize(count);

          for (int i = 0; i < count; i++) {
            double value = D["vec"][i];
            settings.testCase.D[i] = value;
          }
        }
      }
      catch(const SettingNotFoundException &nfex)
      {
        cerr << "No 'D' (diameter) setting in 'dataSet' section in configuration file." << endl;
      }

      // Get Number of Data to generate
      try
      {
        int nData = dataSet.lookup("nData");
        settings.testCase.nData = nData;
      }
      catch(const SettingNotFoundException &nfex)
      {
        cerr << "No 'nData' setting in 'dataSet' section in configuration file." << endl;
      }
    } else if (settings.type == "load") {
      try
      {
        string datafile = dataSet.lookup("file");
        settings.dataInFile = datafile;
      }
      catch(const SettingNotFoundException &nfex)
      {
        cerr << "No 'file' specified for the 'load' setting in 'dataSet' section in configuration file." << endl;
      }

      try
      {
        bool printData = dataSet.lookup("printData");
        settings.printData = printData;
      }
      catch(const SettingNotFoundException &nfex)
      {
        // Ignore
      }
    }
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No 'dataSet' setting in configuration file." << endl;
  }

  // Get the Data Analysis information -----------------------------------------
  try
  {
    const Setting &dataAnalysis = root["dataAnalysis"];

    // Approach selection
    try
    {
      string approach = dataAnalysis.lookup("approach");
      settings.approach = approach;
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'approach' setting in 'dataAnalysis' section in configuration file." << endl;
    }

    // Model selection
    try
    {
      string model = dataAnalysis.lookup("model");
      settings.model = model;
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'model' setting in 'dataAnalysis' section in configuration file." << endl;
    }

    // Get shape parameter list
    try
    {
      const Setting &phiList = dataAnalysis["phi"];
      int nPhi = phiList.getLength();
      settings.phi.resize(nPhi);

      for (int i = 0; i < nPhi; i++) {
        const Setting &phi = phiList[i];

        // Get Name
        try
        {
          string name = phi.lookup("name");
          settings.phi[i].name = name;
        }
        catch(const SettingNotFoundException &nfex)
        {
          string name = "Phi " + to_string(i);
          settings.phi[i].name = name;
        }

        // Get Name written in latex format (optional)
        try
        {
          string latex = phi.lookup("latex");
          settings.phi[i].latex = latex;
        }
        catch(const SettingNotFoundException &nfex) // Default
        {
          settings.phi[i].latex = settings.phi[i].name;
        }

        // Get Number of modes
        try
        {
          int K = phi.lookup("K");
          settings.phi[i].K = K;
        }
        catch(const SettingNotFoundException &nfex) // Default
        {
          settings.phi[i].K = 1;
        }

        // Get parameter vector
        try
        {
          const Setting &vec = phi["vec"];

          bool generate;
          if (vec.lookupValue("generate", generate)
              && generate == true) {

            double first = vec["vec"][0];
            double step  = vec["vec"][1];
            double last  = vec["vec"][2];

            for (double value = first; value <= last; value += step) {
            settings.phi[i].vec.push_back(value);
            }
          }
          else
          {
            int count = vec["vec"].getLength();
            settings.phi[i].vec.resize(count);

            for (int j = 0; j < count; j++) {
              double value = vec["vec"][j];
              settings.phi[i].vec[j] = value;
            }
          }
        }
        catch(const SettingNotFoundException &nfex)
        {
          cerr << "No 'vec' argument in 'phi' setting in 'dataAnalysis' section in configuration file." << endl;
        }

      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      cerr << "No 'phi' setting in 'dataAnalysis' section in configuration file." << endl;
    }

    // Get Mixing Lenght Coefficients (optional)
    try
    {
      const Setting &pi = dataAnalysis["pi"];

      bool generate;
      if (pi.lookupValue("generate", generate)
          && generate == true) {

        double first = pi["vec"][0];
        double step  = pi["vec"][1];
        double last  = pi["vec"][2];

        for (double value = first; value <= last; value += step) {
          settings.pi.push_back(value);
        }

      }
      else
      {
        int count = pi["vec"].getLength();
        settings.pi.resize(count);

        for (int i = 0; i < count; i++) {
          double value = pi["vec"][i];
          settings.pi[i] = value;
        }
      }
    }
    catch(const SettingNotFoundException &nfex)
    {
      // Ignore
    }
  }
  catch(const SettingNotFoundException &nfex)
  {
    cerr << "No 'dataAnalysis' setting in configuration file." << endl;
  }

  return settings;
}

void IO::printProgress(int i, int &prog) {
  if (i == 0 && prog == 0) {
    cout << "0%";
    prog = 25;
  } else if (i > 25 && prog == 25) {
    cout << "----------25%";
    prog = 50;
  } else if (i > 50 && prog == 50) {
    cout << "----------50%";
    prog = 75;
  } else if (i > 75 && prog == 75) {
    cout << "----------75%";
    prog = 100;
  } else if (i == 100) {
    cout << "----------100%" << endl;
  }
}


void IO::printVec(const vector<int>& v) {
  cout << '{';
  for (vector<int>::size_type i = 0; i < v.size()-1; i++) {
    cout << v[i] << ", ";
  }
  cout << v.back() << "}\n";
}


void IO::printVec(const vector<double>& v) {
  cout << '{';
  for (vector<double>::size_type i = 0; i < v.size()-1; i++) {
    cout << v[i] << ", ";
  }
  cout << v.back() << "}\n";
}


void IO::printVec(const vector<float>& v) {
  cout << '{';
  for (vector<float>::size_type i = 0; i < v.size()-1; i++) {
    cout << v[i] << ", ";
  }
  cout << v.back() << "}\n";
}


void IO::printPosterior(const vector<double>& v, const string& file) {
  ofstream out;
  out.open(file);

  for (auto i : v) {
    out << i << '\n';
  }

  out.close();
}
