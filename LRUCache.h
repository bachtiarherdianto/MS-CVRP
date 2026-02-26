/*
 * First developed by: Maria SOTO
 * Applied for CVRP by: Flavien LUCAS
 * Enhanced and upgraded by: Bachtiar HERDIANTO
 */

#ifndef LRUCACHE_H
#define LRUCACHE_H

class LRUCache {  // Least Recently Used Cache
public:
  LRUCache() = default;

  LRUCache(int capacity_, int node_num_)
      : dictSize(capacity_),
        dict(node_num_),
        counter(0) {}

  LRUCache(const LRUCache &other)
      : dictSize(other.dictSize),
        dict(other.dict),
        counter(other.counter),
        right(other.right),
        left(other.left) {}

  LRUCache &operator=(const LRUCache &other) {
    this->dictSize = other.dictSize;
    this->dict = other.dict;
    this->counter = other.counter;
    this->right = other.right;
    this->left = other.left;
    return *this;
  }

  inline void put(const int node) {
    if (node != my::dummy && node != my::depot) {

      if (this->dict[node].used) {
        this->remove(node);
        this->place_top(node);

      } else {  // reached the cache size, evict the least recently used entry

        if (this->counter == this->dictSize) {
          this->remove(this->left);
        } else {
          this->counter++;
        }

        this->place_top(node);  // move the recently accessed entry to top
      }

    }

  }

  void reset() {
    this->counter = 0u;
    int curr = this->right;
    while (curr != my::dummy) {
      const int next = this->dict[curr].next;
      this->dict[curr].used = false;
      this->dict[curr].next = my::dummy;
      this->dict[curr].prev = my::dummy;
      curr = next;
    }

    this->right = my::dummy;
    this->left = my::dummy;
  }

  [[nodiscard]] inline int size() const {
    return this->counter;
  }

  [[nodiscard]] inline int begin() const {
    return this->right;
  }

  [[nodiscard]] inline int get_next(int node) const {
    return this->dict[node].next;
  }

protected:
  class Node {
  public:
    int prev = my::dummy;
    int next = my::dummy;
    bool used = false;
  };

private:
  int dictSize;
  std::vector<Node> dict;
  int counter = 0u;
  int right = my::dummy;
  int left = my::dummy;

  inline void remove(const int node) {
    assert(node != my::dummy);
    assert(this->dict[node].used);
    const int prevEntry = this->dict[node].prev;
    const int nextEntry = this->dict[node].next;

    if (prevEntry == my::dummy) {  // head remove
      this->right = nextEntry;
    } else {
      this->dict[prevEntry].next = nextEntry;
    }

    if (nextEntry == my::dummy) {  // tail remove
      this->left = prevEntry;
    } else {
      this->dict[nextEntry].prev = prevEntry;
    }

    this->dict[node].used = false;
    this->dict[node].prev = my::dummy;
    this->dict[node].next = my::dummy;
  }

  inline void place_top(const int node) {
    assert(node != my::dummy);
    assert(!this->dict[node].used);

    this->dict[node].used = true;
    this->dict[node].next = this->right;

    if (this->right != my::dummy) {
      this->dict[this->right].prev = node;
    }

    this->right = node;
    this->dict[node].prev = my::dummy;

    if (this->left == my::dummy) {
      this->left = this->right;
    }

  }

};

#endif //LRUCACHE_H