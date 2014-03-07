#ifdef HAS_SQLiteDB

#include "../include/SQLiteDB.h"
#include <stdexcept>

namespace skylens {
  SQLiteDB::SQLiteResult::SQLiteResult() : counter(0) { }

  unsigned int SQLiteDB::SQLiteResult::getRowCount() {
    return nrows;
  }
  unsigned int SQLiteDB::SQLiteResult::getFieldCount() {
    return ncols;
  }
  char* SQLiteDB::SQLiteResult::getFieldName(unsigned int i) {
    if (i < ncols)
      return table[i];
    else
      return NULL;
  }
  char** SQLiteDB::SQLiteResult::getRow() {
    if (counter+ncols < (nrows+1)*ncols) {
      return (table + (counter+=ncols));
    } else
      return NULL;
  }
  SQLiteDB::SQLiteResult::~SQLiteResult() {
    if (table != NULL)
      sqlite3_free_table(table);
  }


  SQLiteDB::SQLiteDB() :
    db(NULL), zErrMsg(NULL) {
  }

  void SQLiteDB::connect(std::string dbfile) {
    sqlite3_close(db);
    checkRC(sqlite3_open(dbfile.c_str(), &db));
  }

  SQLiteDB::~SQLiteDB() {
    checkRC(sqlite3_close(db));
  }

  void SQLiteDB::exec(std::string query, int (*callback)(void*,int,char**,char**)) {
    checkRC(sqlite3_exec(db, query.c_str(), callback, 0, &zErrMsg));
  }
  
  SQLiteDB::SQLiteResult SQLiteDB::query(std::string query) {
    SQLiteResult result;
    checkRC(sqlite3_get_table(db, query.c_str(), &result.table, &result.nrows, &result.ncols, &zErrMsg));
    return result;
  }

  void SQLiteDB::checkRC(int rc) {
    if (rc) {
      throw std::runtime_error("SQLiteDB: " +std::string(sqlite3_errmsg(db)));
      sqlite3_close(db);
    }
  }
  
  int SQLiteDB::getDBType() {
    return 0;
  }
}

#endif // HAS_SQLiteDB
