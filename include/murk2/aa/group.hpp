#pragma once

#include <murk2/common/bigint.hpp>
#include <murk2/common/simulation.hpp>

#include <c3lt.hpp>

#include <optional>

namespace murk2::aa {
  struct missing_structure : std::runtime_error {
    missing_structure(const char* msg) : std::runtime_error{msg} {};
  };

  template<template<typename> typename NewStructure, typename OldStructure>
  NewStructure<typename OldStructure::elem_t> const& require_structure(OldStructure const& t) {
    auto ptr = dynamic_cast<NewStructure<typename OldStructure::elem_t> const*>(std::addressof(t));
    if (!ptr)
      throw missing_structure{"The requested operation required extra structure that didn't exist"};
    return *ptr;
  }

  template<typename Element>
  struct monoid;

  template<typename Element>
  struct group;

  template<typename Element>
  struct magma {
    using elem_t = Element;

    virtual Element canonicalise(Element&& elem) const { return elem; }

    virtual Element op(Element const&, Element const&) const = 0;
    virtual void op_mut(Element& a, Element const& b) const {
      a = op(a, b);
    }
    virtual void op_iter_mut(Element& a, bigint const& reps = 2) const {
      if (reps == 0) {
        a = require_structure<monoid>(*this).identity();
        return;
      }

      if (reps < 0) {
        if (!try_invert_mut(a))
          throw missing_structure{"Tried to invert non-invertible element"};

        for (bigint i = -2; i > reps; --i)
          op_mut(a, a);
      }
      else {
        for (bigint i = 2; i < reps; ++i)
          op_mut(a, a);
      }
    }
    virtual Element op_iter(Element const& a, bigint const& reps = 2) const {
      Element ret = a;
      op_iter_mut(ret, reps);
      return ret;
    }

    virtual bool is_invertible(Element const&) const { return false; }
    virtual std::optional<Element> try_invert(Element const&) const { return std::nullopt; }
    virtual std::optional<Element> try_invert(Element&& x) const { return try_invert(x); }
    [[nodiscard]]
    virtual bool try_invert_mut(Element& x) const {
      if (auto res = try_invert(x)) {
        x = std::move(*res);
        return true;
      }
      return false;
    }

    virtual ~magma() = default;
  };

  template<typename Element>
  struct finite_monoid : virtual monoid<Element> {
    // TODO: impl https://www.ams.org/journals/mcom/1998-67-224/S0025-5718-98-00968-5/S0025-5718-98-00968-5.pdf
    virtual bigint order() const = 0;
    virtual bigint element_order(Element const& elem) const {
      bigint ret = 1;
      for (Element acc = elem; !this->is_identity(elem); ++ret, this->op_mut(acc, elem));
      return ret;
    }
  };

  template<typename Element>
  struct monoid : virtual magma<Element> {
    virtual Element identity() const = 0;
    virtual bool is_identity(Element const& a) const { return a == identity(); }

    virtual void op_iter_mut(Element& a, bigint const& reps = 2) const override {
      // General purpose fast multiplication algorithm using doubling

      bigint reps_left = reps;
      if (reps < 0) {
        if (!this->try_invert_mut(a))
          throw missing_structure{"Tried to invert non-invertible element"};
        reps_left = -reps;
      }

      // Simple optimisation
      if (reps_left == 0) {
        a = identity();
        return;
      }
      else if (reps_left == 1)
        return;
      else if (reps_left == 2)
        return this->op_mut(a, a);

      Element cursor = a;
      a = identity();
      for (; reps_left; reps_left >>= 1, this->op_mut(cursor, cursor)) {
        if (mpz_odd_p(reps_left.get_mpz_t()))
          this->op_mut(a, cursor);
      }
    }

    virtual ~monoid() = default;
  };

  template<typename Element>
  struct group;

  template<typename Group, typename Element = typename Group::elem_t>
  class group_element;

  template<typename Element>
  struct group : virtual monoid<Element> {
    inline group_element<group<Element>> element(Element elem) const { return group_element{c3lt::managed(this), std::move(elem)}; }
    inline group_element<group<Element>> operator()(Element elem) const { return element(std::move(elem)); }

    virtual Element invert(Element const& a) const = 0;
    virtual void invert_mut(Element& a) const { a = invert(a); }
    virtual bool is_inverse(Element const& a, Element const& b) const { return this->is_identity(this->op(a, b)); }

    bool is_invertible(Element const&) const override final { return true; }
    std::optional<Element> try_invert(Element const& a) const override final { return invert(a); }
    bool try_invert_mut(Element& a) const override final { invert_mut(a); return true; }

    virtual ~group() = default;
  };

  template<typename Element>
  struct cyclic_group : virtual group<Element> {
    virtual Element generator() const = 0;

    virtual ~cyclic_group() = default;
  };

  template<typename Group, typename Element>
  class group_element {
    static_assert(std::derived_from<Group, group<Element>>);
  public:
    c3lt::safe_ptr <const Group> context;
    Element elem;

  public:
    inline group_element operator-() const& { return context->invert(elem); }
    inline group_element&& operator-() && { context->invert_mut(elem); return std::move(*this); }

    inline group_element operator+(group_element const& b) const& { return {context, context->op(elem, b.elem)}; }
    inline group_element&& operator+(group_element const& b)&& { context->op_mut(elem, b.elem); return std::move(*this); }
    inline group_element& operator+=(group_element const& b) { context->op_mut(elem, b.elem); return *this; }

    inline group_element operator-(group_element const& b) const& { return operator+(-b); }
    inline group_element&& operator-(group_element const& b) && { return operator+(-b); }
    inline group_element& operator-=(group_element const& b) { return operator+=(-b); }

    inline group_element operator*(bigint n) const& { return {context, context->op_iter(elem, n)}; }
    inline group_element&& operator*(bigint n)&& { context->op_iter_mut(elem, n); return std::move(*this); }
    inline group_element& operator*=(bigint n) { context->op_iter_mut(elem, n); return *this; }

  public:
    template<typename Group2>
    inline bool operator==(group_element<Group2, Element> const& other) const { return elem == other.elem; }

  public:
    inline group_element(c3lt::safe_ptr <const Group> context_, Element elem_) : context{context_}, elem{context->canonicalise(std::move(elem_))} {}
  };

  template<typename Group, typename T>
  group_element(c3lt::safe_ptr <Group>, T) -> group_element<Group, typename Group::elem_t>;

  template<typename Group, typename Element>
  inline std::ostream& operator<<(std::ostream& os, group_element<Group, Element> const& x) { return os << x.elem; }

  template<typename Element>
  struct ring;

  template<typename Ring, typename Element = typename Ring::elem_t>
  class ring_element;

  template<typename Element>
  struct ring {
    using elem_t = Element;

    inline ring_element<ring<Element>> element(Element elem) const { return ring_element{c3lt::managed(this), std::move(elem)}; }
    inline ring_element<ring<Element>> operator()(Element elem) const { return element(std::move(elem)); }

    virtual c3lt::safe_ptr <const group<Element>> add() const noexcept = 0;
    virtual c3lt::safe_ptr <const monoid<Element>> ring_mul() const noexcept = 0;

    virtual bool is_zero_divisor(elem_t& x);

    virtual ~ring() = default;
  };

  template<typename Element>
  struct finite_ring : virtual ring<Element> {
    virtual bigint order() const { return require_structure<finite_monoid>(*this->add()).order(); }
  };

  template<typename Ring, typename Element>
  class ring_element {
    static_assert(std::derived_from<Ring, ring<Element>>);
  public:
    c3lt::safe_ptr <const Ring> context;
    Element elem;

  public:
    inline ring_element operator-() const& { return {context, context->add()->invert(elem)}; }
    inline ring_element operator-() && { context->add()->invert_mut(elem); return std::move(*this); }

    inline ring_element operator+(Element const& b) const& { return {context, context->add()->op(elem, b)}; }
    inline ring_element operator+(Element const& b)&& { context->add()->op_mut(elem, b); return std::move(*this); }
    inline ring_element& operator+=(Element const& b) { context->add()->op_mut(elem, b); return *this; }
    inline ring_element operator+(ring_element const& b) const& { return operator+(b.elem); }
    inline ring_element operator+(ring_element const& b)&& { return operator+(b.elem); }
    inline ring_element& operator+=(ring_element const& b) { return operator+(b.elem); }

    inline ring_element operator-(Element const& b) const& { return operator+(context->add()->invert(b)); }
    inline ring_element operator-(Element const& b) && { return operator+(context->add()->invert(b)); }
    inline ring_element& operator-=(Element const& b) { return operator+=(context->add()->invert(b)); }
    inline ring_element operator-(ring_element const& b) const& { return operator-(b.elem); }
    inline ring_element operator-(ring_element const& b) && { return operator-(b.elem); }
    inline ring_element& operator-=(ring_element const& b) {  return operator-=(b.elem); }

    inline ring_element operator*(Element const& b) const& { return {context, context->ring_mul()->op(elem, b)}; }
    inline ring_element operator*(Element const& b)&& { context->ring_mul()->op_mut(elem, b); return std::move(*this); }
    inline ring_element& operator*=(Element const& b) { context->ring_mul()->op_mut(elem, b); return *this; }
    inline ring_element operator*(ring_element const& b) const& { return operator*(b.elem); }
    inline ring_element operator*(ring_element const& b)&& {return operator*(b.elem); }
    inline ring_element& operator*=(ring_element const& b) { return operator*=(b.elem); }

    inline ring_element operator/(ring_element const& b) const& { return operator*(b.invert()); }
    inline ring_element operator/(ring_element const& b)&& { return operator*(b.invert()); }
    inline ring_element& operator/=(ring_element const& b) { return operator*=(b.invert()); }

    inline ring_element operator&(bigint const& n) const& { return {context, context->add()->op_iter(elem, n)}; }
    inline ring_element operator&(bigint const& n)&& { context->add()->op_iter_mut(elem, n); return std::move(*this); }
    inline ring_element& operator&=(bigint const& n) { context->add()->op_iter_mut(elem, n); return *this; }

    inline ring_element operator^(bigint const& n) const& { return {context, context->ring_mul()->op_iter(elem, n)}; }
    inline ring_element operator^(bigint const& n)&& { context->ring_mul()->op_iter_mut(elem, n); return std::move(*this); }
    inline ring_element& operator^=(bigint const& n) { context->ring_mul()->op_iter_mut(elem, n); return *this; }

    inline ring_element invert() const& {
      if (auto x = context->ring_mul()->try_invert(elem))
        return ring_element{context, *x};
      else
        throw missing_structure{"Tried to invert non-invertible element"};
    }
    inline ring_element invert() && {
      if (auto x = context->ring_mul()->try_invert(std::move(elem)))
        return ring_element{context, *x};
      else
        throw missing_structure{"Tried to invert non-invertible element"};
    }
    inline std::optional<ring_element> try_invert() const& {
      if (auto x = context->ring_mul()->try_invert(elem))
        return ring_element{context, std::move(*x)};
      else
        return std::nullopt;
    }
    inline std::optional<ring_element> try_invert() && {
      if (auto x = context->ring_mul()->try_invert(std::move(elem)))
        return ring_element{context, std::move(*x)};
      else
        return std::nullopt;
    }
    [[nodiscard]]
    inline bool try_invert_mut() const& {  return context->try_invert_mut(elem); }

  public:
    template<typename Ring2>
    inline bool operator==(ring_element<Ring2, Element> const& other) const { return elem == other.elem; }

  public:
    inline ring_element(c3lt::safe_ptr <const Ring> context_, Element elem_) : context{context_}, elem{/*context->add()->canonicalise(*/std::move(elem_)/*)*/} {}
  };

  template<typename Ring, typename T>
  ring_element(c3lt::safe_ptr <Ring>, T) -> ring_element<Ring, typename Ring::elem_t>;

  template<typename Field, typename Element = typename Field::elem_t>
  using field_element = ring_element<Field, Element>;

  template<typename Ring, typename Element>
  inline std::ostream& operator<<(std::ostream& os, ring_element<Ring, Element> const& x) { return os << x.elem; }

  template<typename Element>
  class field;

  template<typename Element>
  class field : public virtual ring<Element> {
  public:
    using elem_t = Element;
  private:
    class mul_monoid final : public monoid<Element> {
    private:
      // TODO: check if optimisation worth stupid thread check over init order
      c3lt::safe_ptr <const field<Element>> base;

    public:
      Element op(Element const& a, Element const& b) const {
        if (base->add()->is_identity(a) || base->add()->is_identity(b))
          return base->add()->identity();
        else
          return base->mul()->op(a, b);
      }
      virtual void op_mut(Element& a, Element const& b) const {
        if (base->add()->is_identity(a) || base->add()->is_identity(b))
          a = base->add()->identity();
        else
          return base->mul()->op_mut(a, b);
      }
      virtual void op_iter_mut(Element& a, bigint const& reps = 2) const {
        if (reps < 0 && base->add()->is_identity(a))
          throw missing_structure{"Tried to invert absorbing element"};
        return base->mul()->op_iter_mut(a, std::move(reps));
      }
      virtual Element op_iter(Element const& a, bigint const& reps = 2) const {
        if (reps < 0 && base->add()->is_identity(a))
          throw missing_structure{"Tried to invert absorbing element"};
        return base->mul()->op_iter(a, reps);
      }
      virtual Element op_iter(Element&& a, bigint const& reps = 2) const {
        if (reps < 0 && base->add()->is_identity(a))
          throw missing_structure{"Tried to invert absorbing element"};
        return base->mul()->op_iter(a, std::move(reps));
      }

      virtual bool is_invertible(Element const& a) const { return !base->add()->is_identity(a); }
      virtual std::optional<Element> try_invert(Element const& a) const { return base->add()->is_identity(a) ? std::nullopt : std::optional{base->mul()->invert(a)}; }
      virtual std::optional<Element> try_invert(Element&& a) const { return base->add()->is_identity(a) ? std::nullopt : std::optional{base->mul()->invert(std::move(a))}; }
      virtual bool try_invert_mut(Element& a) const {
        if (base->add()->is_identity(a))
          return false;
        base->mul()->invert_mut(a);
        return true;
      }
      virtual Element identity() const { return base->mul()->identity(); }

    public:
      mul_monoid(c3lt::safe_ptr <const field> f) : base{f} {}
    };

  private:
    mul_monoid mul_;

  public:
    virtual c3lt::safe_ptr <const group<Element>> add() const noexcept override = 0;
    virtual c3lt::safe_ptr <const group<Element>> mul() const noexcept = 0;

    // TODO: work out how to allow copy/move as well
    c3lt::safe_ptr<const monoid<Element>> ring_mul() const noexcept override final { return c3lt::managed(&mul_); }

    field() : mul_{c3lt::managed(this)} {}
    virtual ~field() = default;

    // We keep a ref to mul_, so not possible
    field(field const&) = delete;
    field(field&&) = delete;
  };

  template<typename Element>
  struct finite_field : public virtual field<Element>, public virtual finite_ring<Element> {
    virtual bigint order_prime() const = 0;
    virtual bigint order_exponent() const = 0;
    virtual bigint order() const override {
      bigint ret;
      mpz_pow_ui(ret, order_prime(), order_exponent().get_ui());
      return ret;
    }
    virtual ~finite_field() = default;
  };

  template<typename Monoid>
  class fake_group : public virtual group<typename Monoid::elem_t> {
  private:
    c3lt::safe_ptr <const Monoid> parent_monoid;

  public:
    c3lt::safe_ptr <const Monoid> get_parent_group() const noexcept { return parent_monoid; }

  public:
    virtual typename Monoid::elem_t op(typename Monoid::elem_t const& a, typename Monoid::elem_t const& b) const override { return parent_monoid->op(a, b); }
    virtual void op_mut(typename Monoid::elem_t& a, typename Monoid::elem_t const& b) const override { return parent_monoid->op_mut(a, b); }
    virtual void op_iter_mut(typename Monoid::elem_t& a, bigint const& reps = 2) const override{ return parent_monoid->op_iter_mut(a, reps); }
    virtual typename Monoid::elem_t op_iter(typename Monoid::elem_t const& a, bigint const& reps = 2) const override { return parent_monoid->op_iter(a, reps); }

    virtual typename Monoid::elem_t identity() const override { return parent_monoid->identity(); }

    virtual typename Monoid::elem_t invert(typename Monoid::elem_t const& a) const override {
      if (auto ret = parent_monoid->try_invert(a))
        return ret;
      else
        throw missing_structure{"Fake group called out"};
    }
    virtual void invert_mut(typename Monoid::elem_t& a) const override {
      if (!parent_monoid->try_invert_mut(a))
        throw missing_structure{"Fake group called out"};
    }
    virtual bool is_inverse(typename Monoid::elem_t const& a, typename Monoid::elem_t const& b) const override { return parent_monoid->is_inverse(a, b); }

  public:
    fake_group(c3lt::safe_ptr <const Monoid> parent) : parent_monoid{parent} {}
  };

  template<typename Ring>
  class fake_field : public virtual field<Ring> {
  private:
    c3lt::safe_ptr <const Ring> base;
    fake_group<group<typename Ring::elem_t>> fake_mul;

  public:
    c3lt::safe_ptr <const group<typename Ring::elem_t>> add() const noexcept override final { return base->add(); }
    c3lt::safe_ptr <const group<typename Ring::elem_t>> mul() const noexcept override final { return c3lt::safe_ptr {&fake_mul}; }

  public:
    fake_field(c3lt::safe_ptr <const Ring> base_ring) : base{base_ring}, fake_mul{base_ring->ring_mul()} {};
  };

  /// Used to construct a subgroup from an arbitrary parent group
  template<typename Group>
  class subgroup_helper : public virtual group<typename Group::elem_t> {
  private:
    c3lt::safe_ptr <const Group> parent_group;

  public:
    c3lt::safe_ptr <const Group> get_parent_group() const noexcept { return parent_group; }

  public:
    virtual typename Group::elem_t op(typename Group::elem_t const& a, typename Group::elem_t const& b) const override { return parent_group->op(a, b); }
    virtual void op_mut(typename Group::elem_t& a, typename Group::elem_t const& b) const override { return parent_group->op_mut(a, b); }
    virtual void op_iter_mut(typename Group::elem_t& a, bigint const& reps = 2) const override{ return parent_group->op_iter_mut(a, reps); }
    virtual typename Group::elem_t op_iter(typename Group::elem_t const& a, bigint const& reps = 2) const override { return parent_group->op_iter(a, reps); }

    virtual typename Group::elem_t identity() const override { return parent_group->identity(); }

    virtual typename Group::elem_t invert(typename Group::elem_t const& a) const override { return parent_group->invert(a); }
    virtual void invert_mut(typename Group::elem_t& a) const override { return parent_group->invert_mut(a); }
    virtual bool is_inverse(typename Group::elem_t const& a, typename Group::elem_t const& b) const override { return parent_group->is_inverse(a, b); }

  public:
    subgroup_helper(c3lt::safe_ptr <const Group> parent) : parent_group{parent} {}
  };

  /// Represents the group generated by an element
  template<typename Group, bool specialise_finite = false>
  class cyclic_subgroup;
  template<typename Group>
  class cyclic_subgroup<Group, false> : public subgroup_helper<Group>, public virtual cyclic_group<typename Group::elem_t> {
  protected:
    typename Group::elem_t gen;
    std::optional<typename Group::elem_t> order;

  public:
    typename Group::elem_t generator() const override final { return gen; }

    inline group_element<cyclic_subgroup> element(group_element<Group> const& elem) const { return {c3lt::managed(this), elem.elem}; }
    inline group_element<cyclic_subgroup> element(group_element<Group>&& elem) const { return {c3lt::managed(this), std::move(elem.elem)}; }
    inline group_element<cyclic_subgroup> operator()(group_element<Group> const& elem) const { return {c3lt::managed(this), elem.elem}; }
    inline group_element<cyclic_subgroup> operator()(group_element<Group>&& elem) const { return {c3lt::managed(this), std::move(elem.elem)}; }

  public:
    cyclic_subgroup(c3lt::safe_ptr <const Group> parent, typename Group::elem_t generator_element) : subgroup_helper<Group>{parent}, gen{std::move(generator_element)} {}
    cyclic_subgroup(group_element<Group> generator_element) : subgroup_helper<Group>{generator_element.context}, gen{std::move(generator_element.elem)} {}
  };
  template<typename Group>
  class cyclic_subgroup<Group, true> : public virtual cyclic_subgroup<Group, false>, public virtual finite_monoid<typename Group::elem_t> {
  private:
    typename Group::elem_t order_val;

  public:
    typename Group::elem_t order() const override final { return order_val; }

    inline group_element<cyclic_subgroup> element(group_element<Group> const& elem) const { return {c3lt::managed(this), elem.elem}; }
    inline group_element<cyclic_subgroup> element(group_element<Group>&& elem) const { return {c3lt::managed(this), std::move(elem.elem)}; }
    inline group_element<cyclic_subgroup> operator()(group_element<Group> const& elem) const { return {c3lt::managed(this), elem.elem}; }
    inline group_element<cyclic_subgroup> operator()(group_element<Group>&& elem) const { return {c3lt::managed(this), std::move(elem.elem)}; }

  public:
    cyclic_subgroup(c3lt::safe_ptr <const Group> parent, typename Group::elem_t generator_element) : cyclic_subgroup<Group, false>{parent, std::move(generator_element)}, order_val{parent->element_order(this->gen)} {}
    cyclic_subgroup(c3lt::safe_ptr <const Group> parent, typename Group::elem_t generator_element, bigint order) : cyclic_subgroup<Group, false>{parent, std::move(generator_element)}, order_val{std::move(order)} {}
    cyclic_subgroup(group_element<Group> generator_element) : cyclic_subgroup<Group, false>{std::move(generator_element)}, order_val{this->get_parent_group()->element_order(generator_element)} {}
    cyclic_subgroup(group_element<Group> generator_element, bigint order) : cyclic_subgroup<Group, false>{std::move(generator_element)}, order_val{std::move(order)} {}
  };

  template<typename Group, typename T>
  cyclic_subgroup(c3lt::safe_ptr <Group>, T) -> cyclic_subgroup<Group, false>;
  template<typename Group>
  cyclic_subgroup(group_element<Group>) -> cyclic_subgroup<Group, false>;
}

namespace murk2 {
  template<typename Ring>
  uint64_t simulation_hash(aa::ring_element<Ring> const& elem) {
    return murk2::simulation_hash(elem.elem) & murk2::simulation_hash(elem.context.get());
  }
  template<typename Group>
  uint64_t simulation_hash(aa::group_element<Group> const& elem) {
    return murk2::simulation_hash(elem.elem) & murk2::simulation_hash(elem.context.get());
  }
}
