#pragma once

class settings_manager{
  public:
    settings_manager();

    settings_manager(const settings_manager& s);

    settings_manager(const char* saving_path_);
    
    ~settings_manager();

    void load();

    void load(const char* saving_path_);

    void save() const;

    void save(const char* saving_path_);

    void print() const;

    void set(const char* name, double value);

    void set(const char* name, const char* value);

    double get_double(const char* option) const;

    std::string get_string(const char* option) const;

    void user_input();

    settings_manager& operator=(const settings_manager& s);

  private:
    int number_of_settings;
    std::string saving_path;
    double* values_d;
    std::string* values_s;
    std::string* options;
};
