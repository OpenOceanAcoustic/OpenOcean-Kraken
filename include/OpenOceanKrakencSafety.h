#ifndef OPEN_OCEAN_KRAKENC_SAFETY_H
#define OPEN_OCEAN_KRAKENC_SAFETY_H

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace OpenOceanKrakenc
{
enum class LoadErrorCode
{
    None,
    MissingFile,
    IoError,
    ParseError,
    SchemaError,
    InvalidField,
    UnsupportedCapability
};

struct LoadError
{
    LoadErrorCode code = LoadErrorCode::None;
    std::string path;
    std::string message;
};

struct LoadResult
{
    bool ok = false;
    LoadError error;

    explicit operator bool() const noexcept { return ok; }

    static LoadResult success()
    {
        LoadResult result;
        result.ok = true;
        return result;
    }

    static LoadResult failure(LoadErrorCode code, std::string path,
                              std::string message)
    {
        LoadResult result;
        result.error.code = code;
        result.error.path = std::move(path);
        result.error.message = std::move(message);
        return result;
    }
};

namespace detail
{
struct ResultViewState
{
    std::atomic<std::uint64_t> generation{0};
};
}

template <typename T>
class ArrayView
{
    static_assert(std::is_const_v<T>, "result ArrayView must be read-only");

public:
    using value_type = std::remove_const_t<T>;
    using pointer = T *;
    using reference = T &;
    using const_iterator = pointer;
    using size_type = std::size_t;

    ArrayView() = default;

    ArrayView(pointer data, size_type size,
              std::weak_ptr<const detail::ResultViewState> state,
              std::uint64_t generation)
        : data_(data), size_(size), state_(std::move(state)),
          generation_(generation) {}

    size_type size() const
    {
        validate();
        return size_;
    }
    bool empty() const
    {
        validate();
        return size_ == 0;
    }
    pointer data() const
    {
        validate();
        return data_;
    }
    reference operator[](size_type index) const
    {
        validate();
        return data_[index];
    }
    reference at(size_type index) const
    {
        validate();
        if (index >= size_)
        {
            throw std::out_of_range("OpenOcean-Krakenc result view index is out of range");
        }
        return data_[index];
    }
    const_iterator begin() const
    {
        validate();
        return data_;
    }
    const_iterator end() const
    {
        validate();
        return data_ + size_;
    }

private:
    void validate() const
    {
        const auto state = state_.lock();
        if (!state || state->generation.load(std::memory_order_acquire) != generation_)
        {
            throw std::logic_error("OpenOcean-Krakenc result view was invalidated");
        }
    }

    pointer data_ = nullptr;
    size_type size_ = 0;
    std::weak_ptr<const detail::ResultViewState> state_;
    std::uint64_t generation_ = 0;
};
}

#endif
