extern "C" {
  int get_n_components();

  void get_component_amount(int *, double *);
  void get_component_name(int *, char *);

  void init();
  void minimize();
}
