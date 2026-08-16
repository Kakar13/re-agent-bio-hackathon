"""Pluggable registries for independently versioned scientific providers."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Protocol, runtime_checkable

from re_agent.immuno.contracts import (
    MHCProviderResult,
    MHCRequest,
    ResponseModelRequest,
    ResponseModelResult,
)


@runtime_checkable
class ResponseModelAdapter(Protocol):
    adapter_id: str
    version: str

    def predict(self, request: ResponseModelRequest) -> ResponseModelResult: ...


@runtime_checkable
class MHCProviderAdapter(Protocol):
    provider_id: str
    version: str

    def predict(self, request: MHCRequest) -> MHCProviderResult: ...


class DuplicateProviderError(ValueError):
    pass


class UnknownProviderError(KeyError):
    pass


class ResponseModelRegistry:
    def __init__(self, adapters: Iterable[ResponseModelAdapter] = ()) -> None:
        self._adapters: dict[str, ResponseModelAdapter] = {}
        for adapter in adapters:
            self.register(adapter)

    def register(self, adapter: ResponseModelAdapter) -> None:
        if not isinstance(adapter, ResponseModelAdapter):
            raise TypeError("adapter does not satisfy ResponseModelAdapter")
        if adapter.adapter_id in self._adapters:
            raise DuplicateProviderError(adapter.adapter_id)
        self._adapters[adapter.adapter_id] = adapter

    def get(self, adapter_id: str) -> ResponseModelAdapter:
        try:
            return self._adapters[adapter_id]
        except KeyError as exc:
            raise UnknownProviderError(adapter_id) from exc

    def predict_all(
        self, request: ResponseModelRequest, adapter_ids: Iterable[str] | None = None
    ) -> list[ResponseModelResult]:
        ids = list(adapter_ids) if adapter_ids is not None else sorted(self._adapters)
        return [self.get(adapter_id).predict(request) for adapter_id in ids]

    @property
    def adapter_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._adapters))


class MHCProviderRegistry:
    def __init__(self, providers: Iterable[MHCProviderAdapter] = ()) -> None:
        self._providers: dict[str, MHCProviderAdapter] = {}
        for provider in providers:
            self.register(provider)

    def register(self, provider: MHCProviderAdapter) -> None:
        if not isinstance(provider, MHCProviderAdapter):
            raise TypeError("provider does not satisfy MHCProviderAdapter")
        if provider.provider_id in self._providers:
            raise DuplicateProviderError(provider.provider_id)
        self._providers[provider.provider_id] = provider

    def get(self, provider_id: str) -> MHCProviderAdapter:
        try:
            return self._providers[provider_id]
        except KeyError as exc:
            raise UnknownProviderError(provider_id) from exc

    def predict_all(
        self, request: MHCRequest, provider_ids: Iterable[str] | None = None
    ) -> list[MHCProviderResult]:
        ids = list(provider_ids) if provider_ids is not None else sorted(self._providers)
        return [self.get(provider_id).predict(request) for provider_id in ids]

    @property
    def provider_ids(self) -> tuple[str, ...]:
        return tuple(sorted(self._providers))
