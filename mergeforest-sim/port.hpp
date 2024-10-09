#ifndef MERGEFOREST_SIM_PORT_HPP
#define MERGEFOREST_SIM_PORT_HPP

#include <cstddef>
#include <cstdint>
#include <climits>
#include <cassert>

namespace mergeforest_sim {

using Address = uint64_t;

inline constexpr unsigned mem_transaction_size = 32;
inline constexpr unsigned element_size = 12;
inline constexpr unsigned block_size = 8;
inline constexpr std::size_t block_size_bytes = element_size * block_size;
inline constexpr Address invalid_address = UINT64_MAX;

template<typename Send, typename Recv>
class Port {
  friend Port<Recv, Send>;
public:
  void connect(Port<Recv, Send>* port) {
    other = port;
    port->other = this;
  }

  void transfer() {
    assert(other != nullptr);
    if (!msg_send_valid) return;
    if (other->msg_recv_valid) return;
    other->msg_recv = msg_send;
    other->msg_recv_valid = true;
    msg_send_valid = false;
  }

  bool has_msg_send() const {
    return msg_send_valid;
  }

  bool add_msg_send(const Send& msg) {
    if (msg_send_valid) return false;
    msg_send = msg;
    msg_send_valid = true;
    return true;
  }

  bool msg_received_valid() const {
    return msg_recv_valid;
  }

  Recv get_msg_received() const {
    return msg_recv;
  }

  void clear_msg_received() {
    msg_recv_valid = false;
    msg_recv = {};
  }

  void reset() {
    msg_send = {};
    msg_recv = {};
  }
private:
  Send msg_send {};
  Recv msg_recv {};
  bool msg_send_valid {false};
  bool msg_recv_valid {false};
  Port<Recv, Send>* other {nullptr};
};

struct Empty_Msg {};

struct Mem_Request {
  bool valid() const {
    return address != invalid_address;
  }

  Address address {invalid_address};
  unsigned id {0};
  bool is_write {false};
};

struct Mem_Response {
  bool valid() const {
    return address != invalid_address;
  }

  Address address {invalid_address};
  unsigned id {};
};

namespace mergeforest {
  
struct Prefetched_Row {
  uint32_t B_row_ptr {UINT32_MAX};
  unsigned row_head_ptr {UINT_MAX};
};

struct Cache_Read {
  bool valid() const { return row_ptr != UINT_MAX; }

  unsigned row_ptr {UINT_MAX};
  unsigned id {};
};

struct Cache_Write {
  enum Type {
    write,
    write_last,
    invalid
  };

  bool valid() const { return type != invalid; }
  
  Type type {invalid};
  unsigned num_elements {};
};

struct Cache_Response {
  bool valid() const { return num_elements > 0; }
  unsigned row_ptr {UINT_MAX};
  unsigned num_elements {};
  unsigned id {};
};

} // namespace mergeforest

} // namespace mergeforest_sim

#endif // MERGEFOREST_SIM_PORT_HPP
