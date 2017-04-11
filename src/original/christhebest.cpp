template<typename T1, typename T2, typename T3>
bool myGetline(std::ifstream& stream_, T1& Var1, T2& Var2, T3& Var3) {
  std::string strLine;
  std::istringstream tokenLine;
  if(std::getline(stream_,strLine)) {
    tokenLine.str(strLine);
    tokenLine >> Var1 >> Var2 >> Var3;
    return true;
  }
  else { return false; }
}
